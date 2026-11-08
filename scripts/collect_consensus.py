#!/usr/bin/env python3
"""Per-amplicon consensus.fasta and cluster_sizes.tsv from the polished seeds of one sample.

Inputs: seeds.tsv (consensus, amplicon, cluster_size), the polished FASTA (seed ids
unchanged; a seed missing from it keeps its unpolished sequence from seeds.fasta
next to seeds.tsv) and reads_per_consensus.tsv (consensus, primary mapped reads).
"""

import os


def read_fasta(path):
    seqs, name, parts = {}, None, []
    if not path or not os.path.exists(path):
        return seqs
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def read_table(path):
    """[(consensus, amplicon, cluster_size)] of seeds.tsv."""
    rows = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            if line.strip():
                cid, amp, size = line.rstrip("\n").split("\t")[:3]
                rows.append((cid, amp, int(size)))
    return rows


def read_counts(path):
    counts = {}
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    cid, n = line.rstrip("\n").split("\t")[:2]
                    counts[cid] = int(n)
    return counts


def collect(seeds, polished, fallback, counts, amplicons, cons_paths, sizes_paths):
    by_amp = {amp: [] for amp in amplicons}
    for cid, amp, size in seeds:
        by_amp.setdefault(amp, []).append((cid, size))
    for amp, cons_path, sizes_path in zip(amplicons, cons_paths, sizes_paths):
        os.makedirs(os.path.dirname(cons_path), exist_ok=True)
        with open(cons_path, "w") as out, open(sizes_path, "w") as sz:
            sz.write("consensus\tcluster_size\treads_mapped\n")
            for cid, size in by_amp.get(amp, []):
                seq = polished.get(cid, fallback.get(cid, ""))
                out.write(f">{cid}\n")
                for j in range(0, len(seq), 80):
                    out.write(seq[j:j + 80] + "\n")
                sz.write(f"{cid}\t{size}\t{counts.get(cid, 0)}\n")


def main(smk):
    seeds = read_table(smk.input.table)
    polished = read_fasta(smk.input.polished)
    fallback = read_fasta(os.path.join(os.path.dirname(smk.input.table), "seeds.fasta"))
    counts = read_counts(smk.input.counts)
    collect(seeds, polished, fallback, counts, list(smk.params.amplicons), list(smk.output.cons), list(smk.output.sizes))
    with open(smk.log[0], "w") as log:
        missing = [cid for cid, _, _ in seeds if cid not in polished]
        log.write(f"{len(seeds)} seed(s), {len(polished)} polished, {sum(counts.values())} reads counted\n")
        if missing:
            log.write("not in the polished FASTA (unpolished seed kept): " + ", ".join(missing) + "\n")


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
