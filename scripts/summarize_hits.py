#!/usr/bin/env python3
import os
import sys
import pandas as pd

try:
    hits_file   = snakemake.input.hits
    output_tsv  = snakemake.output.tsv
    min_pid     = float(snakemake.params.min_pid)
    min_qcov    = float(snakemake.params.min_qcov)
    top_delta   = float(snakemake.params.top_delta)
except NameError:
    if len(sys.argv) != 6:
        print("Usage: summarize_hits.py <hits.tsv> <out.tsv> <min_pid> <min_qcov> <top_delta>", file=sys.stderr)
        sys.exit(2)
    hits_file, output_tsv, min_pid, min_qcov, top_delta = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])

def write_empty(path: str):
    pd.DataFrame(columns=["consensus","best_species","pident","qcov","status"]).to_csv(path, sep="\t", index=False)

if not os.path.exists(hits_file) or os.path.getsize(hits_file) == 0:
    write_empty(output_tsv)
    sys.exit(0)

df = pd.read_csv(hits_file, sep="\t", header=None)

if df.empty:
    write_empty(output_tsv)
    sys.exit(0)

if df.shape[1] >= 12:
    cols = ["qseqid","sseqid","pident","alnlen","mismatch","gapopen",
            "qstart","qend","sstart","send","evalue","bits"]
    df = df.iloc[:, :len(cols)].copy()
    df.columns = cols
    has_qcov = False
elif df.shape[1] == 4:
    df.columns = ["qseqid","sseqid","pident","qcov"]
    has_qcov = True
else:
    newcols = ["qseqid","sseqid","pident"] + [f"c{i}" for i in range(df.shape[1]-3)]
    df.columns = newcols
    has_qcov = "qcov" in df.columns

for c in ("pident","qcov","bits","alnlen"):
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

if not has_qcov:
    if {"qstart","qend"}.issubset(df.columns):
        qlen_est = df.groupby("qseqid")["qend"].transform("max").clip(lower=1)
        df["qcov"] = (df["qend"] - df["qstart"] + 1) / qlen_est * 100.0
        has_qcov = True
    else:
        df["qcov"] = 100.0
        has_qcov = True

def parse_species(s):
    if pd.isna(s):
        return "UNK"
    s = str(s)
    return s.split("|species=")[-1] if "|species=" in s else "UNK"

df["species"] = df["sseqid"].map(parse_species)

pid_cut = min_pid * 100.0
qcov_cut = min_qcov * 100.0
df_f = df[(df["pident"] >= pid_cut) & (df["qcov"] >= qcov_cut)].copy()

if df_f.empty:
    write_empty(output_tsv)
    sys.exit(0)

out_rows = []
sort_keys = [("pident", False), ("qcov", False)]
if "bits" in df_f.columns:
    sort_keys.append(("bits", False))
if "alnlen" in df_f.columns:
    sort_keys.append(("alnlen", False))
by = [k for k, _ in sort_keys]
asc = [a for _, a in sort_keys]

for q, sub in df_f.groupby("qseqid"):
    sub = sub.sort_values(by=by, ascending=asc, kind="mergesort").reset_index(drop=True)
    if len(sub) == 0:
        out_rows.append((q, "NO_HIT", None, None, "NO_HIT"))
        continue

    top = sub.iloc[0]
    status = "OK"
    if len(sub) > 1:
        p1 = float(top["pident"])
        p2 = float(sub.iloc[1]["pident"])
        if p1 > 0.0 and (p1 - p2) / p1 < top_delta:
            status = "AMBIGUOUS"

    out_rows.append((
        str(q),
        str(top["species"]),
        float(top["pident"]) / 100.0,
        float(top["qcov"])   / 100.0,
        status
    ))

out = pd.DataFrame(out_rows, columns=["consensus","best_species","pident","qcov","status"])
out.to_csv(output_tsv, sep="\t", index=False)
