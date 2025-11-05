import argparse
import random
import gzip
import logging
import sys


def open_maybe_gzip(path, mode):
    """Open plain or gzipped FASTQ depending on extension."""
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def subsample_fastq(input_fastq, output_fastq, num_reads, seed=None):
    """Subsample a given number of reads from a FASTQ or FASTQ.GZ file."""

    if seed is not None:
        random.seed(seed)
        logging.info(f"Random seed set to {seed}")

    logging.info(f"Loading reads from {input_fastq}...")
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

    logging.info(f"Loaded {len(reads)} reads from input file")

    if num_reads > len(reads):
        raise ValueError(f"Requested {num_reads} reads but input contains only {len(reads)} reads.")

    logging.info(f"Sampling {num_reads} reads...")
    sampled_reads = random.sample(reads, num_reads)

    logging.info(f"Writing {num_reads} sampled reads to {output_fastq}...")
    with open_maybe_gzip(output_fastq, "wt") as f:
        for read in sampled_reads:
            f.writelines(read)

    logging.info(f"Successfully wrote {num_reads} reads to {output_fastq}")


def main():
    parser = argparse.ArgumentParser(description="Subsample reads from a FASTQ or FASTQ.GZ file.")
    parser.add_argument("--input", required=True, help="Input FASTQ(.gz) file")
    parser.add_argument("--output", required=True, help="Output FASTQ(.gz) file")
    parser.add_argument("--num_reads", type=int, required=True, help="Number of reads to sample")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging (DEBUG level)")
    parser.add_argument("--log-file", help="Write logs to file instead of stderr")

    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'

    if args.log_file:
        logging.basicConfig(
            level=log_level,
            format=log_format,
            datefmt=date_format,
            filename=args.log_file,
            filemode='w'
        )
        # Also log to console
        console = logging.StreamHandler(sys.stderr)
        console.setLevel(logging.WARNING)
        console.setFormatter(logging.Formatter(log_format, datefmt=date_format))
        logging.getLogger('').addHandler(console)
    else:
        logging.basicConfig(
            level=log_level,
            format=log_format,
            datefmt=date_format,
            stream=sys.stderr
        )

    try:
        subsample_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            num_reads=args.num_reads,
            seed=args.seed
        )
        return 0
    except ValueError as e:
        logging.error(f"Validation error: {e}")
        return 1
    except (IOError, OSError) as e:
        logging.error(f"File I/O error: {e}")
        return 1
    except Exception as e:
        logging.error(f"Unexpected error: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
