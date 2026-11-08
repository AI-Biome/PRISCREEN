#!/usr/bin/env python3
"""Polish the cluster centroids of one sample into consensus sequences.

All amplicons of a sample are polished in one go, so Medaka (whose start-up
costs tens of seconds) runs once per sample instead of once per amplicon.
Seeds are the centroids with at least ``cluster.min_size`` reads (vsearch
``;size=N``), renamed ``<sample>_<amplicon>_c<i>``. The reads of all bins are
mapped to all seeds with primary alignments only, so every read polishes one
seed, in ``racon_rounds`` rounds of Racon; Medaka then polishes the Racon draft
(unless ``consensus.medaka`` is false). The reads mapped to each consensus in
the final alignment are counted and written next to the cluster size to
``cluster_sizes.tsv`` (reads of clusters below ``min_size`` are thereby
re-assigned to the nearest consensus). Amplicons without seeds get empty files.
"""

import gzip
import os
import re
import shutil
import subprocess
from collections import Counter

SIZE_RE = re.compile(r";size=(\d+)")


def run(cmd, log):
    with open(log, "a") as lf:
        lf.write(f"$ {cmd}\n")
        lf.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=lf, stderr=subprocess.STDOUT, executable="/bin/bash")


def read_fasta(path):
    """[(id, size, seq)]; size from a vsearch ;size= annotation (1 when absent)."""
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


def count_mapped(sam_or_bam):
    """Primary mapped reads per reference from a SAM/BAM file."""
    out = subprocess.run(f"samtools view -F 0x904 {sam_or_bam} | cut -f3", shell=True, check=True,
                         capture_output=True, text=True, executable="/bin/bash").stdout
    return Counter(out.split())


def write_outputs(seeds, polished, counts, cons_paths, sizes_paths, amplicons):
    for amp, cons_path, sizes_path in zip(amplicons, cons_paths, sizes_paths):
        os.makedirs(os.path.dirname(cons_path), exist_ok=True)
        with open(cons_path, "w") as out, open(sizes_path, "w") as sz:
            sz.write("consensus\tcluster_size\treads_mapped\n")
            for cid, size, seq in seeds.get(amp, []):
                s = polished.get(cid, seq)
                out.write(f">{cid}\n")
                for j in range(0, len(s), 80):
                    out.write(s[j:j + 80] + "\n")
                sz.write(f"{cid}\t{size}\t{counts.get(cid, 0)}\n")


def main(smk):
    sample, amplicons = smk.wildcards.sample, list(smk.params.amplicons)
    fqs, cents = list(smk.input.fqs), list(smk.input.cents)
    cons_paths, sizes_paths, log = list(smk.output.cons), list(smk.output.sizes), smk.log[0]
    threads, rounds = int(smk.threads), int(smk.params.rounds)
    use_medaka, model, bacteria = bool(smk.params.medaka), smk.params.model, bool(smk.params.bacteria)
    workdir = smk.params.workdir
    os.makedirs(workdir, exist_ok=True)
    os.makedirs(os.path.dirname(log), exist_ok=True)
    open(log, "w").close()

    seeds = {}
    for amp, cent_path in zip(amplicons, cents):
        recs = read_fasta(cent_path)
        seeds[amp] = [(f"{sample}_{amp}_c{i}", size, seq) for i, (_, size, seq) in enumerate(recs, start=1)]
    n_seeds = sum(len(v) for v in seeds.values())
    if n_seeds == 0:
        write_outputs(seeds, {}, {}, cons_paths, sizes_paths, amplicons)
        with open(log, "a") as lf:
            lf.write("no cluster centroids in any amplicon bin; nothing to polish\n")
        return

    seed_fa = os.path.join(workdir, "seeds.fasta")
    with open(seed_fa, "w") as out:
        for amp in amplicons:
            for cid, _, seq in seeds[amp]:
                out.write(f">{cid}\n{seq}\n")

    # Reads of all non-empty bins in one (multi-member) gzip file.
    reads = os.path.join(workdir, "reads.fq.gz")
    with open(reads, "wb") as out:
        for fq in fqs:
            if not fastq_is_empty(fq):
                with open(fq, "rb") as fh:
                    shutil.copyfileobj(fh, out)

    draft, sam = seed_fa, os.path.join(workdir, "reads_to_draft.sam")
    for r in range(rounds):
        run(f"minimap2 -t {threads} -ax map-ont --secondary=no {draft} {reads} > {sam}", log)
        racon_out = os.path.join(workdir, f"racon{r + 1}.fasta")
        run(f"racon -t {threads} -u {reads} {sam} {draft} > {racon_out}", log)
        draft = racon_out

    if use_medaka:
        medaka_dir = os.path.join(workdir, "medaka")
        shutil.rmtree(medaka_dir, ignore_errors=True)
        opts = f"-m {model}" if model else ("--bacteria" if bacteria else "")
        run(f"medaka_consensus -i {reads} -d {draft} -o {medaka_dir} -t {threads} {opts}", log)
        counts = count_mapped(os.path.join(medaka_dir, "calls_to_draft.bam"))
        draft = os.path.join(medaka_dir, "consensus.fasta")
    else:
        run(f"minimap2 -t {threads} -ax map-ont --secondary=no {draft} {reads} > {sam}", log)
        counts = count_mapped(sam)

    polished = {cid: seq for cid, _, seq in read_fasta(draft)}
    write_outputs(seeds, polished, counts, cons_paths, sizes_paths, amplicons)

    if not bool(smk.params.keep_intermediate):
        for f in [sam, reads, os.path.join(workdir, "medaka", "consensus_probs.hdf"),
                  os.path.join(workdir, "medaka", "calls_to_draft.bam"),
                  os.path.join(workdir, "medaka", "calls_to_draft.bam.bai")]:
            if os.path.exists(f):
                os.remove(f)
        for f in os.listdir(workdir):
            if f.endswith((".mmi", ".fai")):
                os.remove(os.path.join(workdir, f))


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
