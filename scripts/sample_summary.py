#!/usr/bin/env python3
import pandas as pd
import glob
import os

# Get parameters from Snakemake
summaries = snakemake.input.summaries
output_tsv = snakemake.output.tsv

rows = []
for p in summaries:
    amp = p.split("/")[-2]
    if os.path.exists(p) and os.path.getsize(p) > 0:
        df = pd.read_csv(p, sep="\t")
        df.insert(0, "amplicon", amp)
        rows.append(df)

if rows:
    all_df = pd.concat(rows, ignore_index=True)
else:
    all_df = pd.DataFrame(columns=["amplicon","consensus","best_species","pident","qcov","status"])

all_df.to_csv(output_tsv, sep="\t", index=False)
