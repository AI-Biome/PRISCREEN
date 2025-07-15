import argparse
import random


def subsample_fastq(input_fastq, output_fastq, num_reads, seed=None):
    """
    Subsample a given number of reads from a FASTQ file.

    Args:
        input_fastq (str): Path to the input FASTQ file.
        output_fastq (str): Path to save the subsampled FASTQ file.
        num_reads (int): Number of reads to sample.
        seed (int, optional): Random seed for reproducibility.
    """
    if seed is not None:
        random.seed(seed)

    with open(input_fastq, 'r') as f:
        lines = f.readlines()

    reads = [lines[i:i+4] for i in range(0, len(lines), 4)]
    if num_reads > len(reads):
        raise ValueError("Requested more reads than available in input file.")

    sampled_reads = random.sample(reads, num_reads)

    with open(output_fastq, 'w') as f:
        for read in sampled_reads:
            f.writelines(read)


def main():
    parser = argparse.ArgumentParser(description="Subsample a specified number of reads from a FASTQ file.")
    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output FASTQ file")
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

