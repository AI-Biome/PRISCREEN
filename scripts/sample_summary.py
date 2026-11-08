#!/usr/bin/env python3
"""Concatenate the per-amplicon tables of one sample into one TSV (first column: amplicon)."""

import os

import pandas as pd


def concat_tables(paths, columns):
    frames = []
    for p in paths:
        amp = os.path.basename(os.path.dirname(p))
        if os.path.exists(p) and os.path.getsize(p) > 0:
            df = pd.read_csv(p, sep="\t")
            if not df.empty:
                df.insert(0, "amplicon", amp)
                frames.append(df)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["amplicon"] + list(columns))


def main(smk) -> None:
    columns = list(smk.params.columns)
    concat_tables(list(smk.input), columns).to_csv(smk.output.tsv, sep="\t", index=False)


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
