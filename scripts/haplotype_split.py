#!/usr/bin/env python3
"""Split consensus sequences of mixed clusters into haplotypes (one job per sample).

Species variants of a PRISTINE amplicon are often closer to each other than the
raw read error, so read clustering lumps them into one cluster whose consensus is
the majority species. This step recovers the minority:

1. The reads of all bins of a sample are mapped to the polished consensus
   sequences (minimap2, primary alignments only).
2. For every consensus, candidate sites are collected from the alignment:
   substitutions (base counts per position, base quality >= ``min_base_quality``),
   deletions (reads whose deletion starts at a position) and insertions (reads
   with an insertion after a position). A site is a candidate when its
   alternative allele is seen in at least ``min_reads`` reads and in at least
   ``min_allele_frac`` of the reads. A 1-nt indel inside a homopolymer run of
   ``hp_min_run`` or more bases is the typical ONT error and is a *weak* site:
   it needs ``hp_min_allele_frac`` of the reads and never starts a split on its
   own. ``indels: false`` restricts the step to substitutions.
3. Reads are phased by recursive bisection: the strong candidate with the most
   minority reads seeds two groups; every other candidate whose alleles agree
   with that grouping in at least ``min_linkage`` of the reads is kept (random
   errors fall at the same sites in both groups and are dropped); reads are
   re-assigned by majority agreement with the two allele profiles. A split is
   accepted when both groups have at least ``min_reads`` reads and
   ``min_variants`` linked sites. Groups are split again until nothing splits or
   ``max_haplotypes`` is reached.
4. Each haplotype gets the parent consensus with its own alleles at the linked
   sites (substitution, deletion or the most common inserted sequence of its
   reads), named ``<parent>_h<i>`` (largest first). A consensus that does not
   split keeps its name and sequence, so downstream results are unchanged.

Outputs per amplicon mirror ``results/consensus``: ``consensus.fasta`` and
``cluster_sizes.tsv`` (``reads_mapped`` recounted on this alignment; haplotype
rows carry ``parent``, ``n_variants`` and ``separated_by`` = snp, indel or
snp+indel), plus ``variants.tsv`` with every candidate site, its counts,
homopolymer context and the allele of each haplotype.
"""

import gzip
import os
import shutil
import subprocess
from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

BASES = "ACGT"
CODE = {b: i + 1 for i, b in enumerate(BASES)}  # 0 = no call; 1-4 = A, C, G, T
REF, DEL, INS = 1, 5, 6                          # codes in deletion / insertion columns (REF = no indel)
N_CODES = 6
KIND_ORDER = {"ins": 0, "del": 1, "snp": 2}      # edit order at one position (descending positions)


@dataclass
class Params:
    min_reads: int = 5
    min_allele_frac: float = 0.15
    min_depth: int = 20
    min_base_quality: int = 10
    min_linkage: float = 0.8
    min_variants: int = 1
    max_haplotypes: int = 4
    indels: bool = True
    hp_min_run: int = 4
    hp_min_allele_frac: float = 0.2
    max_iter: int = 10


@dataclass
class Site:
    kind: str            # snp, del, ins
    pos: int             # 0-based reference position (ins: insertion after this position)
    ref: str             # reference base (snp), deleted bases (del) or "-" (ins)
    alt: str             # alternative base (snp), "-" (del) or the most common inserted sequence (ins)
    depth: int
    alt_count: int
    hp_run: int          # homopolymer run length at the site (0 = none)
    weak: bool           # 1-nt indel in a homopolymer run of >= hp_min_run bases


@dataclass
class Haplotype:
    reads: np.ndarray                                        # indices into the read list
    alleles: Dict[int, int] = field(default_factory=dict)    # site index -> allele code (linked sites only)
    n_variants: int = 0

    def separated_by(self, sites: Sequence[Site]) -> str:
        kinds = {"indel" if sites[j].kind != "snp" else "snp" for j in self.alleles}
        return "+".join(sorted(kinds)) if kinds else ""


def run(cmd, log):
    with open(log, "a") as lf:
        lf.write(f"$ {cmd}\n")
        lf.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=lf, stderr=subprocess.STDOUT, executable="/bin/bash")


def read_fasta(path):
    """[(id, seq)] of a FASTA file (empty list for missing or empty files)."""
    records, name, parts = [], None, []
    if not path or not os.path.exists(path):
        return records
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.upper())
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def write_fasta(records, out, width=80):
    for name, seq in records:
        out.write(f">{name}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i:i + width] + "\n")


def fastq_is_empty(path):
    with gzip.open(path, "rb") as fh:
        return not fh.read(1)


# ----------------------------------------------------------------------------
# Candidate sites
# ----------------------------------------------------------------------------
def _primary(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate)


def hp_run(seq: str, pos: int, base: str) -> int:
    """Length of the homopolymer run of ``base`` that contains or touches position ``pos``."""
    if not base or base not in BASES:
        return 0
    i = pos
    if seq[pos] != base:
        if pos + 1 < len(seq) and seq[pos + 1] == base:
            i = pos + 1
        elif pos - 1 >= 0 and seq[pos - 1] == base:
            i = pos - 1
        else:
            return 0
    left = i
    while left - 1 >= 0 and seq[left - 1] == base:
        left -= 1
    right = i
    while right + 1 < len(seq) and seq[right + 1] == base:
        right += 1
    return right - left + 1


def base_counts(bam, contig: str, length: int, p: Params) -> np.ndarray:
    """(4, L) counts of A, C, G, T from primary alignments."""
    cov = bam.count_coverage(contig, 0, length, quality_threshold=p.min_base_quality, read_callback=_primary)
    return np.array(cov, dtype=np.int64)


def indel_counts(bam, contig: str, length: int):
    """Deletion starts per position ({pos: Counter(length)}) and insertions after a position ({pos: Counter(seq)})."""
    dels: Dict[int, Counter] = defaultdict(Counter)
    ins: Dict[int, Counter] = defaultdict(Counter)
    for read in bam.fetch(contig):
        if not _primary(read) or read.cigartuples is None:
            continue
        rpos, qpos, seq = read.reference_start, 0, read.query_sequence.upper()
        for op, n in read.cigartuples:
            if op in (0, 7, 8):             # M, =, X
                rpos += n
                qpos += n
            elif op == 2:                   # D
                if rpos < length:
                    dels[rpos][n] += 1
                rpos += n
            elif op == 1:                   # I
                if rpos > read.reference_start and rpos <= length:
                    ins[rpos - 1][seq[qpos:qpos + n]] += 1
                qpos += n
            elif op == 4:                   # S
                qpos += n
            elif op == 3:                   # N
                rpos += n
    return dels, ins


def candidate_sites(seq: str, counts: np.ndarray, dels, ins, p: Params) -> List[Site]:
    """Candidate substitution, deletion and insertion sites of one consensus, sorted by position."""
    sites: List[Site] = []
    depth = counts.sum(axis=0)
    srt = np.sort(counts, axis=0)
    minor = srt[-2]
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(depth > 0, minor / np.maximum(depth, 1), 0.0)
    for pos in np.flatnonzero((depth >= p.min_depth) & (minor >= p.min_reads) & (frac >= p.min_allele_frac)):
        pos = int(pos)
        order = [BASES[i] for i in np.argsort(counts[:, pos])[::-1]]
        alt = next((b for b in order if b != seq[pos]), order[0])
        sites.append(Site("snp", pos, seq[pos], alt, int(depth[pos]), int(counts[CODE[alt] - 1, pos]), 0, False))
    if p.indels:
        for pos, lengths in dels.items():
            lengths = Counter(lengths)
            n = sum(lengths.values())
            total = int(depth[pos]) + n
            if total < p.min_depth or n < p.min_reads:
                continue
            length = lengths.most_common(1)[0][0]
            run = hp_run(seq, pos, seq[pos])
            weak = length == 1 and run >= p.hp_min_run
            if n / total >= (p.hp_min_allele_frac if weak else p.min_allele_frac):
                sites.append(Site("del", pos, seq[pos:pos + length], "-", total, n, run, weak))
        for pos, seqs in ins.items():
            seqs = Counter(seqs)
            n = sum(seqs.values())
            total = int(depth[pos])
            if total < p.min_depth or n < p.min_reads:
                continue
            alt = seqs.most_common(1)[0][0]
            run = hp_run(seq, pos, alt) if len(alt) == 1 else 0
            weak = len(alt) == 1 and run >= p.hp_min_run
            if n / total >= (p.hp_min_allele_frac if weak else p.min_allele_frac):
                sites.append(Site("ins", pos, "-", alt, total, n, run, weak))
    sites.sort(key=lambda s: (s.pos, KIND_ORDER[s.kind]))
    return sites


def read_alleles(bam, contig: str, sites: Sequence[Site], p: Params):
    """Read names, an (n_reads, n_sites) matrix of allele codes (0 = no call) and the inserted
    sequence per (site index, read index) for insertion sites."""
    snp_j = {s.pos: j for j, s in enumerate(sites) if s.kind == "snp"}
    del_j = {s.pos: j for j, s in enumerate(sites) if s.kind == "del"}
    ins_j = {s.pos: j for j, s in enumerate(sites) if s.kind == "ins"}
    positions = sorted({s.pos for s in sites})
    names, rows = [], []
    ins_seqs: Dict[int, Dict[int, str]] = defaultdict(dict)
    for read in bam.fetch(contig):
        if not _primary(read):
            continue
        ri = len(names)
        names.append(read.query_name)
        row = np.zeros(len(sites), dtype=np.uint8)
        if sites and read.cigartuples is not None:
            seq, quals = read.query_sequence.upper(), read.query_qualities
            rpos, qpos = read.reference_start, 0
            for op, n in read.cigartuples:
                if op in (0, 7, 8):
                    for pos in positions[bisect_left(positions, rpos):bisect_left(positions, rpos + n)]:
                        q = qpos + (pos - rpos)
                        j = snp_j.get(pos)
                        if j is not None and (quals is None or quals[q] >= p.min_base_quality):
                            row[j] = CODE.get(seq[q], 0)
                        j = del_j.get(pos)
                        if j is not None:
                            row[j] = REF
                        j = ins_j.get(pos)
                        if j is not None:
                            row[j] = REF
                    rpos += n
                    qpos += n
                elif op == 2:
                    j = del_j.get(rpos)
                    if j is not None:
                        row[j] = DEL
                    rpos += n
                elif op == 1:
                    j = ins_j.get(rpos - 1) if rpos > read.reference_start else None
                    if j is not None:
                        row[j] = INS
                        ins_seqs[j][ri] = seq[qpos:qpos + n]
                    qpos += n
                elif op == 4:
                    qpos += n
                elif op == 3:
                    rpos += n
        rows.append(row)
    matrix = np.vstack(rows) if rows else np.zeros((0, len(sites)), dtype=np.uint8)
    return names, matrix, ins_seqs


# ----------------------------------------------------------------------------
# Phasing
# ----------------------------------------------------------------------------
def _allele_counts(sub: np.ndarray) -> np.ndarray:
    """(N_CODES, n_sites) counts of the allele codes 1..N_CODES."""
    return np.stack([(sub == c).sum(axis=0) for c in range(1, N_CODES + 1)])


def _majority(sub: np.ndarray) -> np.ndarray:
    """Majority allele code per site (0 where no read has a call)."""
    counts = _allele_counts(sub)
    maj = counts.argmax(axis=0) + 1
    maj[counts.sum(axis=0) == 0] = 0
    return maj.astype(np.uint8)


def bisect(matrix: np.ndarray, idx: np.ndarray, weak: np.ndarray, p: Params):
    """Split the reads idx into two haplotypes. Returns (idx_a, idx_b, linked site indices) or None.

    weak: boolean array per site; weak sites need hp_min_allele_frac and cannot seed a split."""
    sub = matrix[idx]
    if sub.shape[1] == 0 or len(idx) < 2 * p.min_reads:
        return None
    counts = _allele_counts(sub)
    srt = np.sort(counts, axis=0)
    minor, depth = srt[-2], counts.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(depth > 0, minor / np.maximum(depth, 1), 0.0)
    need = np.where(weak, p.hp_min_allele_frac, p.min_allele_frac)
    active = np.flatnonzero((minor >= p.min_reads) & (frac >= need) & (depth >= p.min_depth))
    strong = [j for j in active if not weak[j]]
    if not strong:
        return None
    seed = int(max(strong, key=lambda j: minor[j]))
    order = np.argsort(counts[:, seed])[::-1]
    major_code, minor_code = order[0] + 1, order[1] + 1
    assign = np.full(len(idx), -1, dtype=np.int8)
    assign[sub[:, seed] == major_code] = 0
    assign[sub[:, seed] == minor_code] = 1

    def linked_sites(assign):
        """Candidates whose alleles agree with the current grouping, and the two allele profiles."""
        in_a, in_b = assign == 0, assign == 1
        prof_a, prof_b = _majority(sub[in_a]), _majority(sub[in_b])
        new_linked = []
        for j in active:
            a, b = prof_a[j], prof_b[j]
            if a == 0 or b == 0 or a == b:
                continue
            calls = sub[:, j]
            called = (calls != 0) & (assign >= 0)
            if called.sum() == 0:
                continue
            expected = np.where(assign == 0, a, b)
            concordance = (calls[called] == expected[called]).mean()
            if concordance >= p.min_linkage:
                new_linked.append(int(j))
        return new_linked, prof_a, prof_b

    linked = [seed]
    for _ in range(p.max_iter):
        if (assign == 0).sum() < p.min_reads or (assign == 1).sum() < p.min_reads:
            return None
        linked, prof_a, prof_b = linked_sites(assign)
        if not linked:
            return None
        # re-assign every read by agreement with the two profiles over the linked sites
        calls = sub[:, linked]
        has_call = calls != 0
        score_a = ((calls == prof_a[linked]) & has_call).sum(axis=1)
        score_b = ((calls == prof_b[linked]) & has_call).sum(axis=1)
        new_assign = np.full(len(idx), -1, dtype=np.int8)
        new_assign[score_a > score_b] = 0
        new_assign[score_b > score_a] = 1
        if np.array_equal(new_assign, assign):
            break
        assign = new_assign
    else:
        linked, _, _ = linked_sites(assign)     # max_iter reached: linked sites of the final grouping
        if not linked:
            return None

    in_a, in_b = assign == 0, assign == 1
    n_a, n_b = int(in_a.sum()), int(in_b.sum())
    if min(n_a, n_b) < p.min_reads or len(linked) < p.min_variants or all(weak[j] for j in linked):
        return None
    if min(n_a, n_b) / (n_a + n_b) < p.min_allele_frac:
        return None
    return idx[in_a], idx[in_b], sorted(linked)


def phase(matrix: np.ndarray, p: Params, weak: Optional[np.ndarray] = None) -> List[Haplotype]:
    """Recursively split all reads (rows of matrix) into haplotypes; a single haplotype when nothing splits."""
    if weak is None:
        weak = np.zeros(matrix.shape[1], dtype=bool)
    all_idx = np.arange(matrix.shape[0])
    leaves: List[Haplotype] = [Haplotype(all_idx)]
    open_groups = [0]
    while open_groups and len(leaves) < p.max_haplotypes:
        gi = open_groups.pop(0)
        res = bisect(matrix, leaves[gi].reads, weak, p)
        if res is None:
            continue
        idx_a, idx_b, linked = res
        parent = leaves[gi]
        for k, sub_idx in enumerate((idx_a, idx_b)):
            hap = Haplotype(sub_idx, dict(parent.alleles))
            maj = _majority(matrix[sub_idx])
            for j in linked:
                if maj[j]:
                    hap.alleles[j] = int(maj[j])
            if k == 0:
                leaves[gi] = hap
                open_groups.append(gi)
            else:
                leaves.append(hap)
                open_groups.append(len(leaves) - 1)
    for hap in leaves:
        hap.n_variants = len(hap.alleles)
    if len(leaves) == 1:
        leaves[0].alleles, leaves[0].n_variants = {}, 0
    leaves.sort(key=lambda h: -len(h.reads))
    return leaves


def allele_name(site: Site, code: int) -> str:
    if site.kind == "snp":
        return BASES[code - 1] if 1 <= code <= 4 else "."
    return {REF: "ref", DEL: "del", INS: "ins"}.get(code, ".")


def apply_alleles(seq: str, sites: Sequence[Site], alleles: Dict[int, int], ins_choice: Optional[Dict[int, str]] = None) -> str:
    """Parent consensus with the haplotype's alleles: substitutions, deletions and insertions
    (``ins_choice``: site index -> inserted sequence; default the site's most common one)."""
    edits = sorted(alleles.items(), key=lambda kv: (sites[kv[0]].pos, KIND_ORDER[sites[kv[0]].kind]), reverse=True)
    s = seq
    for j, code in edits:
        site = sites[j]
        if site.kind == "snp" and 1 <= code <= 4:
            s = s[:site.pos] + BASES[code - 1] + s[site.pos + 1:]
        elif site.kind == "del" and code == DEL:
            s = s[:site.pos] + s[site.pos + len(site.ref):]
        elif site.kind == "ins" and code == INS:
            insert = (ins_choice or {}).get(j, site.alt)
            s = s[:site.pos + 1] + insert + s[site.pos + 1:]
    return s


# ----------------------------------------------------------------------------
# Per-sample driver
# ----------------------------------------------------------------------------
def split_consensus(bam, cid: str, seq: str, p: Params):
    """(records, size_rows, variant_rows, reads_unphased) for one consensus sequence.

    records: [(name, sequence)]; size_rows: [(name, reads, parent, n_variants, separated_by)];
    variant_rows: one per candidate site (see write_outputs for the columns)."""
    counts = base_counts(bam, cid, len(seq), p)
    dels, ins = indel_counts(bam, cid, len(seq)) if p.indels else ({}, {})
    sites = candidate_sites(seq, counts, dels, ins, p)
    names, matrix, ins_seqs = read_alleles(bam, cid, sites, p)
    n_reads = len(names)
    weak = np.array([s.weak for s in sites], dtype=bool)
    haps = phase(matrix, p, weak) if sites else [Haplotype(np.arange(n_reads))]

    if len(haps) == 1:
        records = [(cid, seq)]
        sizes = [(cid, n_reads, "", 0, "")]
    else:
        records, sizes = [], []
        for i, hap in enumerate(haps, start=1):
            name = f"{cid}_h{i}"
            choice = {}
            for j, code in hap.alleles.items():
                if sites[j].kind == "ins" and code == INS:
                    seqs_here = Counter(ins_seqs[j][r] for r in hap.reads if r in ins_seqs[j])
                    if seqs_here:
                        choice[j] = seqs_here.most_common(1)[0][0]
            records.append((name, apply_alleles(seq, sites, hap.alleles, choice)))
            sizes.append((name, len(hap.reads), cid, hap.n_variants, hap.separated_by(sites)))
    n_phased = sum(len(h.reads) for h in haps) if len(haps) > 1 else n_reads
    variants = []
    for j, site in enumerate(sites):
        hap_alleles = ";".join(f"h{i}={allele_name(site, h.alleles[j]) if j in h.alleles else '.'}"
                               for i, h in enumerate(haps, start=1)) if len(haps) > 1 else ""
        variants.append((cid, site.pos + 1, site.kind, site.ref, site.alt, site.depth, site.alt_count,
                         site.hp_run, int(site.weak), hap_alleles))
    return records, sizes, variants, n_reads - n_phased


def write_outputs(amplicon, parents, out_cons, out_sizes, out_vars, results):
    os.makedirs(os.path.dirname(out_cons), exist_ok=True)
    with open(out_cons, "w") as fa, open(out_sizes, "w") as sz, open(out_vars, "w") as va:
        sz.write("consensus\tcluster_size\treads_mapped\tparent\tn_variants\tseparated_by\treads_unphased\n")
        va.write("consensus\tposition\ttype\tref\talt\tdepth\talt_count\thomopolymer_run\tweak\thaplotype_alleles\n")
        for cid, cluster_size in parents:
            records, sizes, variants, unphased = results[cid]
            write_fasta(records, fa)
            for k, (name, reads, parent, n_var, sep) in enumerate(sizes):
                sz.write(f"{name}\t{cluster_size}\t{reads}\t{parent}\t{n_var}\t{sep}\t{unphased if (parent and k == 0) else 0}\n")
            for row in variants:
                va.write("\t".join(str(x) for x in row) + "\n")


def read_cluster_sizes(path):
    """[(consensus, cluster_size)] in file order."""
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            out.append((row["consensus"], int(row.get("cluster_size", 0))))
    return out


def main(smk):
    import pysam

    sample, amplicons = smk.wildcards.sample, list(smk.params.amplicons)
    fqs, cons_in, sizes_in = list(smk.input.fqs), list(smk.input.cons), list(smk.input.sizes)
    cons_out, sizes_out, vars_out = list(smk.output.cons), list(smk.output.sizes), list(smk.output.variants)
    log, workdir, threads = smk.log[0], smk.params.workdir, int(smk.threads)
    p = Params(**{k: type(getattr(Params, k))(v) for k, v in dict(smk.params.phasing).items()})
    os.makedirs(workdir, exist_ok=True)
    os.makedirs(os.path.dirname(log), exist_ok=True)
    open(log, "w").close()

    parents = {amp: read_cluster_sizes(s) for amp, s in zip(amplicons, sizes_in)}
    seqs = {}
    for amp, c in zip(amplicons, cons_in):
        for cid, seq in read_fasta(c):
            seqs[cid] = seq
    results = {}
    if seqs:
        ref = os.path.join(workdir, "consensus_all.fasta")
        with open(ref, "w") as out:
            write_fasta(seqs.items(), out)
        reads = os.path.join(workdir, "reads.fq.gz")
        with open(reads, "wb") as out:
            for fq in fqs:
                if not fastq_is_empty(fq):
                    with open(fq, "rb") as fh:
                        shutil.copyfileobj(fh, out)
        bam = os.path.join(workdir, "reads_to_consensus.bam")
        run(f"minimap2 -t {threads} -ax map-ont --secondary=no {ref} {reads} | "
            f"samtools sort -@ {min(4, threads)} -m 1G -o {bam} - && samtools index {bam}", log)
        with pysam.AlignmentFile(bam) as af, open(log, "a") as lf:
            for cid, seq in seqs.items():
                results[cid] = split_consensus(af, cid, seq, p)
                records, sizes, variants, unphased = results[cid]
                sep = ", ".join(f"{name}: {s}" for name, _, parent, _, s in sizes if parent)
                lf.write(f"{cid}: {len(variants)} candidate site(s), {len(records)} haplotype(s)"
                         f"{' (' + sep + ')' if sep else ''}, {unphased} read(s) unphased\n")
        for f in (reads, ref + ".fai"):
            if os.path.exists(f):
                os.remove(f)
        if not bool(smk.params.keep_bam):
            for f in (bam, bam + ".bai"):
                os.remove(f)
    for amp, oc, os_, ov in zip(amplicons, cons_out, sizes_out, vars_out):
        for cid, _ in parents[amp]:
            results.setdefault(cid, ([(cid, seqs.get(cid, ""))], [(cid, 0, "", 0, "")], [], 0))
        write_outputs(amp, parents[amp], oc, os_, ov, results)


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
