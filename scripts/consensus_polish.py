"""
Consensus polishing script for Snakemake
"""
import os
import subprocess
import tempfile
import gzip

def deduplicate_fastq(input_fq, output_fq, amplicon_filter=None):
    """Deduplicate FASTQ reads by sequence ID, keeping first occurrence

    Args:
        input_fq: Input FASTQ file (gzipped)
        output_fq: Output FASTQ file (uncompressed)
        amplicon_filter: If provided, only keep reads matching this amplicon name
    """
    seen_ids = set()
    with gzip.open(input_fq, "rt") as infile, open(output_fq, "w") as outfile:
        while True:
            header = infile.readline()
            if not header:
                break
            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()

            # Extract read ID (remove @ and everything after first space)
            read_id = header[1:].split()[0].strip()

            # Filter by amplicon if specified
            if amplicon_filter:
                # Extract amplicon from read ID (format: sample_amplicon_...)
                parts = read_id.split("_")
                if len(parts) >= 2:
                    read_amplicon = parts[1]  # e.g., "amp1" from "Test1_amp1_A_0001"
                    if read_amplicon != amplicon_filter:
                        continue  # Skip reads from other amplicons

            if read_id not in seen_ids:
                seen_ids.add(read_id)
                outfile.write(header)
                outfile.write(seq)
                outfile.write(plus)
                outfile.write(qual)

def run_racon_rounds(fastq_path, initial_consensus, output_path, rounds, threads):
    """Run multiple rounds of racon polishing"""
    current_consensus = initial_consensus

    for i in range(rounds):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_cons:
            tmp_cons_path = tmp_cons.name
            if i == 0:
                # First round: use initial consensus
                subprocess.check_call(f"cp {current_consensus} {tmp_cons_path}", shell=True)
            else:
                # Subsequent rounds: use previous output
                with open(tmp_cons_path, "w") as f:
                    with open(current_consensus) as src:
                        f.write(src.read())

        # Align reads to consensus
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sam", delete=False) as tmp_sam:
            tmp_sam_path = tmp_sam.name

        subprocess.check_call(
            f"minimap2 -ax map-ont -t {threads} {tmp_cons_path} {fastq_path} > {tmp_sam_path}",
            shell=True
        )

        # Run racon
        if i == rounds - 1:
            # Final round: write to output
            subprocess.check_call(
                f"racon -t {threads} {fastq_path} {tmp_sam_path} {tmp_cons_path} > {output_path}",
                shell=True
            )
        else:
            # Intermediate round: write to temp file
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_out:
                tmp_out_path = tmp_out.name
            subprocess.check_call(
                f"racon -t {threads} {fastq_path} {tmp_sam_path} {tmp_cons_path} > {tmp_out_path}",
                shell=True
            )
            current_consensus = tmp_out_path

        # Cleanup
        os.unlink(tmp_cons_path)
        os.unlink(tmp_sam_path)

def main(snakemake):
    """Generate consensus sequences with racon polishing"""

    # Create output directory
    os.makedirs(
        f"results/consensus/{snakemake.wildcards.sample}/{snakemake.wildcards.amplicon}",
        exist_ok=True
    )

    # Get amplicon name for filtering
    amplicon = snakemake.wildcards.amplicon

    if snakemake.params.clustered:
        # Deduplicate reads first and filter by amplicon
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as tmp_fq:
            tmp_fq_path = tmp_fq.name

        deduplicate_fastq(snakemake.input.fq, tmp_fq_path, amplicon_filter=amplicon)

        # Run racon polishing on cluster centroids
        run_racon_rounds(
            tmp_fq_path,
            snakemake.input.cents,
            snakemake.output.cons,
            snakemake.params.rounds,
            snakemake.threads
        )

        # Cleanup
        os.unlink(tmp_fq_path)
    else:
        # Use first read as initial consensus, then polish with racon
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_cons:
            tmp_cons_path = tmp_cons.name

        with gzip.open(snakemake.input.fq, "rt") as fh:
            h = fh.readline().strip()
            s = fh.readline().strip()
            if not h.startswith("@"):
                raise RuntimeError("FASTQ appears empty or malformed.")

            with open(tmp_cons_path, "w") as out:
                out.write(">consensus\n" + s + "\n")

        # Deduplicate reads and filter by amplicon
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as tmp_fq:
            tmp_fq_path = tmp_fq.name

        deduplicate_fastq(snakemake.input.fq, tmp_fq_path, amplicon_filter=amplicon)

        # Run racon polishing
        run_racon_rounds(
            tmp_fq_path,
            tmp_cons_path,
            snakemake.output.cons,
            snakemake.params.rounds,
            snakemake.threads
        )

        # Cleanup
        os.unlink(tmp_cons_path)
        os.unlink(tmp_fq_path)

if __name__ == "__main__":
    main(snakemake)
