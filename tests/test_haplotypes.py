"""Haplotype splitting: phasing on synthetic allele matrices and on a real minimap2 alignment."""
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pysam
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import haplotype_split as hs  # noqa: E402
from simulate import apply_edits, simulate_reads, substitute, write_fastq  # noqa: E402

P = hs.Params(min_reads=5, min_allele_frac=0.15, min_depth=20, min_base_quality=10, min_linkage=0.8)


def matrix_from_haplotypes(haps, n_per_hap, n_noise_pos, rng, error=0.02, missing=0.05):
    """Allele matrix of reads from several haplotypes (strings over ACGT at the variant positions)
    plus noise positions where every read has the same base except for random errors."""
    rows, truth = [], []
    for h, hap in enumerate(haps):
        for _ in range(n_per_hap[h]):
            row = []
            for base in hap + "A" * n_noise_pos:
                r = rng.random()
                if r < missing:
                    row.append(0)
                elif r < missing + error:
                    row.append(hs.CODE[rng.choice([b for b in "ACGT" if b != base])])
                else:
                    row.append(hs.CODE[base])
            rows.append(row)
            truth.append(h)
    return np.array(rows, dtype=np.uint8), np.array(truth)


def purity(hap_list, truth):
    return [np.bincount(truth[h.reads]).max() / len(h.reads) for h in hap_list]


def test_two_haplotypes_three_snps():
    rng = random.Random(1)
    m, truth = matrix_from_haplotypes(["AAA", "CCC"], [60, 25], n_noise_pos=4, rng=rng)
    haps = hs.phase(m, P)
    assert len(haps) == 2
    assert sorted(len(h.reads) for h in haps) == pytest.approx([25, 60], abs=3)
    assert min(purity(haps, truth)) > 0.95
    assert all(h.n_variants == 3 for h in haps)
    assert set(haps[0].alleles.values()) == {hs.CODE["A"]} and set(haps[1].alleles.values()) == {hs.CODE["C"]}


def test_single_snp_still_splits_with_min_variants_1():
    rng = random.Random(2)
    m, truth = matrix_from_haplotypes(["G", "T"], [50, 20], n_noise_pos=3, rng=rng)
    assert len(hs.phase(m, P)) == 2
    strict = hs.Params(**{**P.__dict__, "min_variants": 2})
    assert len(hs.phase(m, strict)) == 1


def test_three_haplotypes_recursive():
    rng = random.Random(3)
    m, truth = matrix_from_haplotypes(["AAAA", "CCAA", "AACC"], [60, 30, 30], n_noise_pos=3, rng=rng)
    haps = hs.phase(m, P)
    assert len(haps) == 3
    assert min(purity(haps, truth)) > 0.9
    assert len(hs.phase(m, hs.Params(**{**P.__dict__, "max_haplotypes": 2}))) == 2


def test_noise_only_does_not_split():
    rng = random.Random(4)
    m, _ = matrix_from_haplotypes(["AAAAAA"], [120], n_noise_pos=4, rng=rng, error=0.12)
    assert len(hs.phase(m, P)) == 1


def test_candidate_sites_thresholds_and_homopolymers():
    seq = "AAAAACGTAC"
    counts = np.zeros((4, 10), dtype=int)
    counts[0] = [100] * 10                       # A everywhere
    counts[1, 1] = 20                            # C at pos 1: 20/120 passes
    counts[1, 2] = 4                             # 4 < min_reads
    counts[1, 3] = 30
    dels = {0: {1: 15}, 6: {1: 20}, 8: {2: 20}}  # pos 0: 1-nt in AAAAA = weak (15/115 < 0.2); pos 6: strong; pos 8: 2-nt
    ins = {4: {"A": 40}, 7: {"GATC": 30}}        # after pos 4: A into AAAAA = weak, 40/100 passes; GATC strong
    sites = hs.candidate_sites(seq, counts, dels, ins, P)
    assert [(s.kind, s.pos, s.weak) for s in sites] == [
        ("snp", 1, False), ("snp", 3, False), ("ins", 4, True), ("del", 6, False), ("ins", 7, False), ("del", 8, False)]
    assert [s for s in sites if s.kind == "ins" and s.pos == 4][0].hp_run == 5
    assert [s.ref for s in sites if s.kind == "del"] == ["G", "AC"]
    no_indels = hs.Params(**{**P.__dict__, "indels": False})
    assert [s.kind for s in hs.candidate_sites(seq, counts, dels, ins, no_indels)] == ["snp", "snp"]


def test_apply_alleles():
    sites = [hs.Site("snp", 1, "C", "T", 0, 0, 0, False), hs.Site("del", 3, "TA", "-", 0, 0, 0, False),
             hs.Site("ins", 5, "-", "GG", 0, 0, 0, False), hs.Site("snp", 6, "G", "A", 0, 0, 0, False)]
    assert hs.apply_alleles("ACGTACGT", sites, {0: hs.CODE["T"], 1: hs.DEL, 2: hs.INS, 3: hs.CODE["A"]}) == "ATGCGGAT"
    assert hs.apply_alleles("ACGTACGT", sites, {1: hs.REF, 2: hs.INS}, {2: "TTT"}) == "ACGTACTTTGT"


def test_weak_sites_cannot_split_alone():
    rng = random.Random(5)
    m, truth = matrix_from_haplotypes(["G", "T"], [50, 20], n_noise_pos=2, rng=rng)
    weak = np.array([True, False, False])
    assert len(hs.phase(m, P, weak)) == 1
    assert len(hs.phase(m, P, np.zeros(3, dtype=bool))) == 2


@pytest.fixture(scope="module")
def mixed_alignment(tmp_path_factory):
    """Reads of two variants (3 SNPs apart) mapped to the majority consensus; plus a pure second consensus."""
    rng = random.Random(7)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    var = substitute(ref, [100, 300, 500], rng)
    other = "".join(rng.choice("ACGT") for _ in range(500))
    d = tmp_path_factory.mktemp("hap")
    cons = d / "cons.fasta"
    cons.write_text(f">c1\n{ref}\n>c2\n{other}\n")
    reads = simulate_reads({"maj": (ref, 70), "min": (var, 30), "oth": (other, 40)}, rng)
    fq = d / "reads.fq.gz"
    write_fastq(reads, str(fq))
    bam = d / "reads.bam"
    subprocess.run(f"minimap2 -ax map-ont --secondary=no {cons} {fq} 2>/dev/null | samtools sort -o {bam} - && samtools index {bam}",
                   shell=True, check=True)
    return cons, bam, ref, var, other


def test_split_consensus_on_alignment(mixed_alignment):
    cons, bam, ref, var, other = mixed_alignment
    with pysam.AlignmentFile(str(bam)) as af:
        records, sizes, variants, unphased = hs.split_consensus(af, "c1", ref, P)
        rec2, sizes2, var2, unphased2 = hs.split_consensus(af, "c2", other, P)
    assert [r[0] for r in records] == ["c1_h1", "c1_h2"]
    assert {records[0][1], records[1][1]} == {ref, var}
    reads = {name: n for name, n, *_ in sizes}
    assert abs(reads["c1_h1"] - 70) <= 3 and abs(reads["c1_h2"] - 30) <= 3 and unphased <= 3
    assert sorted(v[1] for v in variants if v[9]) == [101, 301, 501]
    assert all(v[2] == "snp" for v in variants if v[9])
    assert sizes[0][4] == "snp"
    assert rec2 == [("c2", other)] and sizes2 == [("c2", 40, "", 0, "")] and unphased2 == 0


def test_main_like_run(mixed_alignment, tmp_path):
    """The outputs written by main(): per-amplicon FASTA, cluster_sizes.tsv and variants.tsv."""
    cons, bam, ref, var, other = mixed_alignment
    with pysam.AlignmentFile(str(bam)) as af:
        results = {cid: hs.split_consensus(af, cid, seq, P) for cid, seq in (("c1", ref), ("c2", other))}
    out_fa, out_sz, out_va = tmp_path / "consensus.fasta", tmp_path / "cluster_sizes.tsv", tmp_path / "variants.tsv"
    hs.write_outputs("ampX", [("c1", 90), ("c2", 40)], str(out_fa), str(out_sz), str(out_va), results)
    names = [l[1:].strip() for l in out_fa.read_text().splitlines() if l.startswith(">")]
    assert names == ["c1_h1", "c1_h2", "c2"]
    rows = [l.split("\t") for l in out_sz.read_text().splitlines()]
    assert rows[0] == ["consensus", "cluster_size", "reads_mapped", "parent", "n_variants", "separated_by", "reads_unphased"]
    assert rows[1][0] == "c1_h1" and rows[1][1] == "90" and rows[1][3] == "c1" and rows[1][4] == "3" and rows[1][5] == "snp"
    assert rows[2][6] == "0"       # reads_unphased is given once per cluster, on the first haplotype row
    assert rows[3] == ["c2", "40", "40", "", "0", "", "0"]
    assert len(out_va.read_text().splitlines()) >= 4


@pytest.fixture(scope="module")
def indel_alignment(tmp_path_factory):
    """Two variants that differ by three indels only: a 1-nt deletion outside homopolymers, a 1-nt insertion
    into a homopolymer run of 5, and a 4-nt insertion; reads with ~2 % error and homopolymer indel errors."""
    rng = random.Random(11)
    ref = list(rng.choice("ACGT") for _ in range(600))
    ref[200:205] = "GGGGG"
    ref[90], ref[91], ref[92] = "A", "C", "T"           # no homopolymer around the deletion
    ref = "".join(ref)
    var = apply_edits(ref, [(91, "del", 1), (202, "ins", "G"), (400, "ins", "GATC")])
    d = tmp_path_factory.mktemp("hapindel")
    cons = d / "cons.fasta"
    cons.write_text(f">c1\n{ref}\n")
    reads = simulate_reads({"maj": (ref, 70), "min": (var, 30)}, rng, sub=0.01, ins=0.005, dele=0.005, hp_indel=0.15)
    fq = d / "reads.fq.gz"
    write_fastq(reads, str(fq))
    bam = d / "reads.bam"
    subprocess.run(f"minimap2 -ax map-ont --secondary=no {cons} {fq} 2>/dev/null | samtools sort -o {bam} - && samtools index {bam}",
                   shell=True, check=True)
    return bam, ref, var


def test_split_by_indels_only(indel_alignment):
    bam, ref, var = indel_alignment
    with pysam.AlignmentFile(str(bam)) as af:
        records, sizes, variants, unphased = hs.split_consensus(af, "c1", ref, P)
    assert [r[0] for r in records] == ["c1_h1", "c1_h2"]
    assert records[0][1] == ref and records[1][1] == var
    reads = {name: n for name, n, _, _, _ in sizes}
    assert abs(reads["c1_h1"] - 70) <= 3 and abs(reads["c1_h2"] - 30) <= 3 and unphased <= 3
    assert sizes[1][4] == "indel"
    linked = [(v[2], v[1]) for v in variants if "h2=del" in v[9] or "h2=ins" in v[9]]
    assert ("del", 92) in linked and ("ins", 401) in linked
    with pysam.AlignmentFile(str(bam)) as af:
        no_indels = hs.Params(**{**P.__dict__, "indels": False})
        assert len(hs.split_consensus(af, "c1", ref, no_indels)[0]) == 1
