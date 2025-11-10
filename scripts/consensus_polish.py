#!/usr/bin/env python3
import os
import subprocess
import tempfile
import gzip

# Get parameters from Snakemake
sample = snakemake.wildcards.sample
amplicon = snakemake.wildcards.amplicon
racon = snakemake.input.racon
fq = snakemake.input.fq
cents = snakemake.input.cents
output_cons = snakemake.output.cons
threads = snakemake.threads
rounds = snakemake.params.rounds
model = snakemake.params.model
clustered = snakemake.params.clustered

os.makedirs(f"results/consensus/{sample}/{amplicon}", exist_ok=True)

# Check if input FASTQ is empty
fq_is_empty = False
try:
    with gzip.open(fq, "rt") as fh:
        first_line = fh.readline()
        if not first_line or not first_line.startswith("@"):
            fq_is_empty = True
except:
    fq_is_empty = True

# If FASTQ is empty, create empty output and exit
if fq_is_empty:
    with open(output_cons, "w") as out:
        pass  # Create empty file
    exit(0)

seed_fa = f"results/consensus/{sample}/{amplicon}/seed.fa"
if clustered:
    # Check if centroids file exists and is not empty
    if not os.path.exists(cents) or os.path.getsize(cents) == 0:
        with open(output_cons, "w") as out:
            pass  # Create empty file
        exit(0)
    subprocess.check_call(f"cp {cents} {seed_fa}", shell=True)
else:
    with gzip.open(fq, "rt") as fh, open(seed_fa, "w") as out:
        h = fh.readline().strip()
        s = fh.readline().strip()
        fh.readline()
        fh.readline()
        if not h.startswith("@"):
            raise RuntimeError("FASTQ appears empty or malformed for seed.")
        out.write(">seed\n" + s + "\n")

tmp_sam = f"results/consensus/{sample}/{amplicon}/aln.sam"

def map_and_sort(template_fa):
    subprocess.check_call(
        f"minimap2 -t {threads} -ax map-ont --secondary=no {template_fa} {fq} > {tmp_sam}",
        shell=True
    )

map_and_sort(seed_fa)
for r in range(rounds):
    racon_out = f"results/consensus/{sample}/{amplicon}/racon{r}.fa"
    subprocess.check_call(
        f"{racon} -t {threads} {fq} {tmp_sam} {seed_fa} > {racon_out}",
        shell=True
    )
    seed_fa = racon_out
    map_and_sort(seed_fa)

medaka_dir = f"results/consensus/{sample}/{amplicon}/medaka"
subprocess.check_call(f"rm -rf {medaka_dir}", shell=True)
subprocess.check_call(
    f"medaka_consensus -i {fq} -d {seed_fa} -o {medaka_dir} -t {threads} -m {model}",
    shell=True
)
subprocess.check_call(f"cp {medaka_dir}/consensus.fasta {output_cons}", shell=True)
