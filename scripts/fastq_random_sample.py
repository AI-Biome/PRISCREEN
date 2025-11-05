import argparse
import random
import gzip


def open_maybe_gzip(path, mode):
    """Open plain or gzipped FASTQ depending on extension."""
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def subsample_fastq(input_fastq, output_fastq, num_reads, seed=None):
    """Subsample a given number of reads from a FASTQ or FASTQ.GZ file."""

    if seed is not None:
        random.seed(seed)

    reads = []
    with open_maybe_gzip(input_fastq, "rt") as f:
        while True:
            r1 = f.readline()
            if not r1:
                break  # EOF
            r2 = f.readline()
            r3 = f.readline()
            r4 = f.readline()
            reads.append([r1, r2, r3, r4])

    if num_reads > len(reads):
        raise ValueError(f"Requested {num_reads} reads but input contains only {len(reads)} reads.")

    sampled_reads = random.sample(reads, num_reads)


    with open_maybe_gzip(output_fastq, "wt") as f:
        for read in sampled_reads:
            f.writelines(read)


def main():
    parser = argparse.ArgumentParser(description="Subsample reads from a FASTQ or FASTQ.GZ file.")
    parser.add_argument("--input", required=True, help="Input FASTQ(.gz) file")
    parser.add_argument("--output", required=True, help="Output FASTQ(.gz) file")
    parser.add_argument("--num_reads", type=int, required=True, help="Number of reads to sample")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed")

    args = parser.parse_args()

    subsample_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        num_reads=args.num_reads,
        seed=args.seed
    )


if __name__ == "__main__":
    main()
