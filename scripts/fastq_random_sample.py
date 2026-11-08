#!/usr/bin/env python3
"""Subsample a fixed number of reads from a FASTQ(.gz) file in one streaming pass.

Reservoir sampling keeps only ``num_reads`` records in memory, so files of any
size can be subsampled. With ``--seed`` the result is reproducible.

    python scripts/fastq_random_sample.py --input in.fastq.gz --output out.fastq.gz --num_reads 10000 --seed 42
"""
import argparse
import gzip
import random


def open_maybe_gzip(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def fastq_records(fh):
    while True:
        header = fh.readline()
        if not header:
            return
        seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"malformed FASTQ record near: {header.strip()!r}")
        yield header, seq, plus, qual


def reservoir_sample(records, k, rng):
    """k records chosen uniformly at random from an iterator (Algorithm R)."""
    sample = []
    for n, rec in enumerate(records):
        if n < k:
            sample.append(rec)
        else:
            j = rng.randrange(n + 1)
            if j < k:
                sample[j] = rec
    return sample, n + 1 if sample else 0


def subsample_fastq(input_fastq, output_fastq, num_reads, seed=None):
    rng = random.Random(seed)
    with open_maybe_gzip(input_fastq, "rt") as fin:
        sample, total = reservoir_sample(fastq_records(fin), num_reads, rng)
    if total < num_reads:
        raise ValueError(f"requested {num_reads} reads but the input has only {total}")
    with open_maybe_gzip(output_fastq, "wt") as fout:
        for rec in sample:
            fout.writelines(rec)
    return total


def main():
    ap = argparse.ArgumentParser(description="Subsample reads from a FASTQ or FASTQ.GZ file (streaming).")
    ap.add_argument("--input", required=True, help="input FASTQ(.gz)")
    ap.add_argument("--output", required=True, help="output FASTQ(.gz)")
    ap.add_argument("--num_reads", type=int, required=True, help="number of reads to keep")
    ap.add_argument("--seed", type=int, default=None, help="random seed")
    args = ap.parse_args()
    total = subsample_fastq(args.input, args.output, args.num_reads, args.seed)
    print(f"kept {args.num_reads} of {total} reads")


if __name__ == "__main__":
    main()
