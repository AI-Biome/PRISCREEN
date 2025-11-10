#!/usr/bin/env python3
import pandas as pd
import os

# Get parameters from Snakemake
hits_file = snakemake.input.hits
output_tsv = snakemake.output.tsv
min_pid = snakemake.params.min_pid
min_qcov = snakemake.params.min_qcov
top_delta = snakemake.params.top_delta

cols = ["qseqid","sseqid","pident","qcov","alnlen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bits"]
df = pd.read_csv(hits_file, sep="\t", header=None, names=cols)

if df.empty:
    out = pd.DataFrame(columns=["consensus","best_species","pident","qcov","status"])
    out.to_csv(output_tsv, sep="\t", index=False)
    exit(0)

def parse_species(s):
    return s.split("|species=")[-1] if "|species=" in s else "UNK"

df["species"] = df["sseqid"].map(parse_species)

df_f = df[(df["pident"] >= min_pid*100) & (df["qcov"] >= min_qcov*100)].copy()

out_rows = []
for q, sub in df_f.groupby("qseqid"):
    sub = sub.sort_values(["pident","qcov","bits"], ascending=[False, False, False]).reset_index(drop=True)
    if len(sub)==0:
        out_rows.append((q, "NO_HIT", None, None, "NO_HIT"))
        continue
    top = sub.iloc[0]
    status = "OK"
    if len(sub) > 1:
        p1 = top["pident"]
        p2 = sub.iloc[1]["pident"]
        if p1 > 0 and (p1 - p2)/p1 < top_delta:
            status = "AMBIGUOUS"
    out_rows.append((q, top["species"], float(top["pident"])/100.0, float(top["qcov"])/100.0, status))

out = pd.DataFrame(out_rows, columns=["consensus","best_species","pident","qcov","status"])
out.to_csv(output_tsv, sep="\t", index=False)
