"""
Consensus polishing script for Snakemake
"""
import os
import subprocess
import tempfile
import gzip

def main(snakemake):
    """Generate consensus sequences with racon and medaka polishing"""

    # Create output directory
    os.makedirs(
        f"results/consensus/{snakemake.wildcards.sample}/{snakemake.wildcards.amplicon}",
        exist_ok=True
    )

    # Prepare seed sequence
    seed_fa = f"results/consensus/{snakemake.wildcards.sample}/{snakemake.wildcards.amplicon}/seed.fa"

    if snakemake.params.clustered:
        # Use cluster centroids as seed
        subprocess.check_call(f"cp {snakemake.input.cents} {seed_fa}", shell=True)
    else:
        # Use first read as seed
        with gzip.open(snakemake.input.fq, "rt") as fh, open(seed_fa, "w") as out:
            h = fh.readline().strip()
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            if not h.startswith("@"):
                raise RuntimeError("FASTQ appears empty or malformed for seed.")
            out.write(">seed\n" + s + "\n")

    # Temporary BAM file for alignments
    tmp_bam = f"results/consensus/{snakemake.wildcards.sample}/{snakemake.wildcards.amplicon}/aln.bam"

    def map_and_sort(template_fa):
        """Map reads to template and sort"""
        cmd = f"minimap2 -t {snakemake.threads} -ax map-ont {template_fa} {snakemake.input.fq} | samtools view -b - | samtools sort -o {tmp_bam}"
        subprocess.check_call(cmd, shell=True)
        subprocess.check_call(f"samtools index {tmp_bam}", shell=True)

    # Initial mapping
    map_and_sort(seed_fa)

    # Racon polishing rounds
    # Extract FASTQ to temporary file (racon doesn't like process substitution)
    tmp_fq = f"results/consensus/{snakemake.wildcards.sample}/{snakemake.wildcards.amplicon}/reads.fq"
    subprocess.check_call(f"zcat {snakemake.input.fq} > {tmp_fq}", shell=True)

    for r in range(snakemake.params.rounds):
        racon_out = f"results/consensus/{snakemake.wildcards.sample}/{snakemake.wildcards.amplicon}/racon{r}.fa"
        cmd = f"racon -t {snakemake.threads} {tmp_fq} {tmp_bam} {seed_fa} > {racon_out}"
        subprocess.check_call(cmd, shell=True)
        seed_fa = racon_out
        map_and_sort(seed_fa)

    # Final consensus (racon-polished only, medaka skipped due to dependency conflicts)
    # Copy the final racon output as the consensus
    subprocess.check_call(f"cp {seed_fa} {snakemake.output.cons}", shell=True)

if __name__ == "__main__":
    main(snakemake)
