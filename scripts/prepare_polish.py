#!/usr/bin/env python3
"""Prepare the polishing of one sample: seeds of all amplicons and one read file.

Seeds are the cluster centroids with at least ``cluster.min_size`` reads (vsearch
``;size=N`` annotation), renamed ``<sample>_<amplicon>_c<i>``. ``seeds.tsv``
keeps the amplicon and cluster size of every seed; ``reads.fq.gz`` is the
concatenation of all non-empty bins (a multi-member gzip file, which minimap2,
Racon and Medaka read). Pure Python, so it runs in any Python container.
"""

import gzip
import os
import re
import shutil

SIZE_RE = re.compile(r";size=(\d+)")


def read_centroids(path):
    """[(name, size, seq)] of a vsearch centroid FASTA (size from ;size=, 1 when absent)."""
    records, name, size, parts = [], None, 1, []
    if not path or not os.path.exists(path):
        return records
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, size, "".join(parts)))
                head = line[1:].split()[0]
                m = SIZE_RE.search(head)
                size = int(m.group(1)) if m else 1
                name, parts = SIZE_RE.sub("", head), []
            else:
                parts.append(line)
    if name is not None:
        records.append((name, size, "".join(parts)))
    return records


def fastq_is_empty(path):
    with gzip.open(path, "rb") as fh:
        return not fh.read(1)


def prepare(sample, amplicons, cent_paths, fq_paths, seeds_fa, seeds_tsv, reads_out):
    n_seeds = 0
    with open(seeds_fa, "w") as fa, open(seeds_tsv, "w") as tsv:
        tsv.write("consensus\tamplicon\tcluster_size\n")
        for amp, path in zip(amplicons, cent_paths):
            for i, (_, size, seq) in enumerate(read_centroids(path), start=1):
                cid = f"{sample}_{amp}_c{i}"
                fa.write(f">{cid}\n{seq}\n")
                tsv.write(f"{cid}\t{amp}\t{size}\n")
                n_seeds += 1
    with open(reads_out, "wb") as out:
        if n_seeds:
            for fq in fq_paths:
                if not fastq_is_empty(fq):
                    with open(fq, "rb") as fh:
                        shutil.copyfileobj(fh, out)
    return n_seeds


def main(smk):
    n = prepare(smk.wildcards.sample, list(smk.params.amplicons), list(smk.input.cents), list(smk.input.fqs),
                smk.output.seeds, smk.output.table, smk.output.reads)
    with open(smk.log[0], "w") as log:
        log.write(f"{n} seed(s)\n" if n else "no cluster centroids in any amplicon bin; nothing to polish\n")


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
