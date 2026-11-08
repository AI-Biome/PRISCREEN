"""Tiny ONT-like read simulator for the tests (substitutions, insertions, deletions; no quality model)."""
import gzip
import random
from typing import Dict, List, Sequence, Tuple

BASES = "ACGT"


def apply_edits(seq: str, edits: Sequence[Tuple[int, str, str]]) -> str:
    """seq with (0-based position, kind, payload) edits: ("del", length), ("ins", sequence inserted after pos),
    ("sub", base). Applied from the right so positions refer to the input sequence."""
    s = seq
    for pos, kind, payload in sorted(edits, key=lambda e: -e[0]):
        if kind == "del":
            s = s[:pos] + s[pos + int(payload):]
        elif kind == "ins":
            s = s[:pos + 1] + str(payload) + s[pos + 1:]
        elif kind == "sub":
            s = s[:pos] + str(payload) + s[pos + 1:]
    return s


def substitute(seq: str, positions: Sequence[int], rng: random.Random) -> str:
    """seq with a different base at every 0-based position."""
    s = list(seq)
    for pos in positions:
        s[pos] = rng.choice([b for b in BASES if b != s[pos]])
    return "".join(s)


def simulate_read(template: str, rng: random.Random, sub: float = 0.01, ins: float = 0.003, dele: float = 0.003,
                  hp_indel: float = 0.0, hp_min_run: int = 4) -> str:
    """ONT-like read: uniform substitutions and indels plus, with probability ``hp_indel`` per homopolymer
    run of at least ``hp_min_run`` bases, one extra base or one base less in that run."""
    out = []
    n = len(template)
    i = 0
    while i < n:
        base = template[i]
        run_end = i
        while run_end + 1 < n and template[run_end + 1] == base:
            run_end += 1
        run_len = run_end - i + 1
        hp_change = 0
        if run_len >= hp_min_run and rng.random() < hp_indel:
            hp_change = rng.choice((-1, 1))
        for k in range(run_len + hp_change):
            b = base
            r = rng.random()
            if r < dele:
                continue
            if r < dele + ins:
                out.append(rng.choice(BASES))
            if rng.random() < sub:
                b = rng.choice([x for x in BASES if x != b])
            out.append(b)
        i = run_end + 1
    return "".join(out)


def simulate_reads(templates: Dict[str, Tuple[str, int]], rng: random.Random, **error) -> List[Tuple[str, str]]:
    """[(name, sequence)] for {prefix: (template, n_reads)}."""
    reads = []
    for prefix, (template, n) in templates.items():
        for i in range(n):
            reads.append((f"{prefix}_{i}", simulate_read(template, rng, **error)))
    return reads


def write_fastq(reads: Sequence[Tuple[str, str]], path: str, mode: str = "wt") -> None:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, mode) as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
