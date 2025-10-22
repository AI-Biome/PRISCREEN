import argparse
import random
import gzip
from typing import Optional

DNA_ALPHABET = ["A", "C", "G", "T", "N"]

def mutate_base(base: str) -> str:
    """Return a random nucleotide different from the original (case‑preserving)."""
    is_lower = base.islower()
    choices = [b for b in DNA_ALPHABET if b.upper() != base.upper()]
    new_base = random.choice(choices)
    return new_base.lower() if is_lower else new_base


def open_maybe_gzip(path: str, mode: str = "rt"):
    """Open plain or gzipped FASTQ depending on extension."""
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def simulate_errors(
    input_fastq: str,
    output_fastq: str,
    error_rate: float,
    seed: Optional[int] = None,
    low_quality_char: str = "!",
    preserve_quality: bool = False,
):
    """Introduce substitution errors into FASTQ reads at the given per‑base error_rate.

    Parameters
    ----------
    input_fastq : str
        Input FASTQ file path (optionally .gz).
    output_fastq : str
        Output FASTQ file path (same compression as input).
    error_rate : float
        Probability of mutating each base (0–1).
    seed : int, optional
        Random seed for reproducibility.
    low_quality_char : str
        ASCII character to assign to mutated bases (ignored if preserve_quality=True).
    preserve_quality : bool
        If True, keep original quality symbols even at mutated positions.
    """
    if not 0 <= error_rate <= 1:
        raise ValueError("error_rate must be between 0 and 1")
    if seed is not None:
        random.seed(seed)

    with open_maybe_gzip(input_fastq, "rt") as fin, open_maybe_gzip(output_fastq, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break  # EOF
            seq = fin.readline().rstrip()
            plus = fin.readline()
            qual = fin.readline().rstrip()

            if not (header.startswith("@") and plus.startswith("+")):
                raise ValueError("Malformed FASTQ record detected around line: " + header)

            seq_list = list(seq)
            qual_list = list(qual)

            # mutate sequence & optionally quality
            for i in range(len(seq_list)):
                if random.random() < error_rate:
                    seq_list[i] = mutate_base(seq_list[i])
                    if not preserve_quality:
                        qual_list[i] = low_quality_char
            mutated_seq = "".join(seq_list)
            mutated_qual = "".join(qual_list)

            fout.write(header)
            fout.write(mutated_seq + "\n")
            fout.write(plus)
            fout.write(mutated_qual + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Introduce a user‑defined substitution error rate into a FASTQ file."\
    )
    parser.add_argument("--input", required=True, help="Input FASTQ (.fastq or .fastq.gz)")
    parser.add_argument("--output", required=True, help="Output FASTQ (.fastq or .fastq.gz)")
    parser.add_argument("--error_rate", type=float, required=True, help="Per‑base error rate, e.g. 0.01 = 1%")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    parser.add_argument(
        "--low_quality_char",
        default="!",
        help="ASCII char to assign to mutated bases (ignored with --preserve_quality)",
    )
    parser.add_argument(
        "--preserve_quality",
        action="store_true",
        help="Keep original quality symbols even where bases are mutated",
    )

    args = parser.parse_args()

    simulate_errors(
        input_fastq=args.input,
        output_fastq=args.output,
        error_rate=args.error_rate,
        seed=args.seed,
        low_quality_char=args.low_quality_char,
        preserve_quality=args.preserve_quality,
    )


if __name__ == "__main__":
    main()

