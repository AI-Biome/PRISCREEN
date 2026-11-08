#!/usr/bin/env python3
"""Per-read species support for panel amplicons from the binning alignment.

Consensus-based identification cannot separate species whose amplicons differ
by only a few sites when both occur in one sample: their reads fall into the
same cluster and polish into one consensus. This script complements it with a
read-level view. For every read with a primary alignment to a panel amplicon,
the alignment score (AS) of its primary hit is compared with the best secondary
hit to another species of the same amplicon (minimap2 -N 10 --secondary=yes):

* ``reads``: reads whose primary alignment is to this species,
* ``reads_unique``: of those, reads that score strictly better than any other
  species (ties, e.g. identical variants, are not unique),
* ``reads_margin_ge_<n>``: reads that score at least ``min_margin`` better.

The counts are written per amplicon and species. They are a vote, not an
identification: a species that gets many unique votes is likely present even
when its consensus is shared with a relative.
"""

import sys
from collections import defaultdict

import pandas as pd


def species_of(ref: str):
    """(amplicon, species) of a mapping-reference name 'amp|species=X|...'; species None for rRNA bins."""
    fields = ref.split("|")
    amp = fields[0]
    for f in fields[1:]:
        if f.startswith("species="):
            return amp, f[len("species="):]
    return amp, None


def parse_sam(lines, min_margin: int, panel_amplicons=None):
    """Count votes from SAM lines (no header) in any order (the BAM is coordinate-sorted, so the
    secondary records of a read may precede its primary). Only amplicons in ``panel_amplicons``
    (all when None) are counted."""
    primary = {}                      # qname -> (amp, species, AS)
    secondaries = defaultdict(list)   # qname -> [(amp, species, AS)]
    for line in lines:
        if not line or line[0] == "@":
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 11:
            continue
        flag, ref = int(cols[1]), cols[2]
        if flag & 0x4 or flag & 0x800 or ref == "*":
            continue
        amp, sp = species_of(ref)
        if sp is None or (panel_amplicons is not None and amp not in panel_amplicons):
            continue
        score = None
        for t in cols[11:]:
            if t.startswith("AS:i:"):
                score = int(t[5:])
                break
        if score is None:
            continue
        q = cols[0]
        if flag & 0x100:
            secondaries[q].append((amp, sp, score))
        else:
            primary[q] = (amp, sp, score)
    best_other = {}                   # qname -> best AS of another species of the same amplicon
    for q, (amp, sp, _) in primary.items():
        others = [sc for a, s2, sc in secondaries.get(q, ()) if a == amp and s2 != sp]
        if others:
            best_other[q] = max(others)
    votes = defaultdict(lambda: [0, 0, 0])
    for q, (amp, sp, score) in primary.items():
        other = best_other.get(q, -1)
        v = votes[(amp, sp)]
        v[0] += 1
        if other < 0 or score > other:
            v[1] += 1
        if other < 0 or score - other >= min_margin:
            v[2] += 1
    rows = [dict(amplicon=a, species=s, reads=v[0], reads_unique=v[1], **{f"reads_margin_ge_{min_margin}": v[2]})
            for (a, s), v in sorted(votes.items())]
    cols = ["amplicon", "species", "reads", "reads_unique", f"reads_margin_ge_{min_margin}"]
    df = pd.DataFrame(rows, columns=cols)
    if not df.empty:
        df["fraction_of_amplicon_reads"] = df["reads"] / df.groupby("amplicon")["reads"].transform("sum")
    else:
        df["fraction_of_amplicon_reads"] = []
    return df


def main(smk) -> None:
    import subprocess
    min_margin = int(smk.params.min_margin)
    proc = subprocess.Popen(["samtools", "view", smk.input.bam], stdout=subprocess.PIPE, text=True)
    df = parse_sam(proc.stdout, min_margin, set(smk.params.panel_amplicons))
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("samtools view failed")
    df.to_csv(smk.output.tsv, sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
