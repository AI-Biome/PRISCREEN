#!/usr/bin/env python3
import os
import pandas as pd

paths = snakemake.input.flags
out_tsv = snakemake.output.tsv

rows = []
for p in paths:
    amp = p.split("/")[-2]
    if os.path.exists(p) and os.path.getsize(p) > 0:
        df = pd.read_csv(p, sep="\t")
        if "consensus" not in df.columns:
            continue
        df.insert(0, "amplicon", amp)
        rows.append(df)

if rows:
    all_df = pd.concat(rows, ignore_index=True)
else:
    all_df = pd.DataFrame(columns=[
        "amplicon","consensus","best_species","pident","qcov","status",
        "read_count","rel_abundance","is_rare","is_novel","flag_reason","blast_url"
    ])

all_df.to_csv(out_tsv, sep="\t", index=False)