#!/usr/bin/env python3
import os
import urllib.parse
import pandas as pd

summary_path = snakemake.input.summary
cons_fa = snakemake.input.cons
uc_path = snakemake.input.uc
out_tsv = snakemake.output.tsv
out_fa = snakemake.output.fasta

min_reads = int(snakemake.params.min_reads)
min_rel = float(snakemake.params.min_rel)
novel_pid = float(snakemake.params.novel_pid)
novel_qcov = float(snakemake.params.novel_qcov)
unclass_is_novel = bool(snakemake.params.unclass_is_novel)

def read_fasta(path):
    seqs = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return seqs
    name = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            seqs[name] = "".join(parts)
    return seqs

def parse_uc_cluster_sizes(path):
    sizes = {}
    if not path or (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return sizes
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] != "H":
                continue
            centroid = cols[9]
            sizes[centroid] = sizes.get(centroid, 0) + 1
    for c in list(sizes.keys()):
        sizes[c] = sizes[c] + 1
    return sizes

def blast_url(seq):
    q = urllib.parse.quote(seq)
    return "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn&QUERY=" + q

if (not os.path.exists(summary_path)) or os.path.getsize(summary_path) == 0:
    pd.DataFrame(columns=[
        "consensus","best_species","pident","qcov","status",
        "read_count","rel_abundance",
        "is_rare","is_novel","flag_reason",
        "blast_url"
    ]).to_csv(out_tsv, sep="\t", index=False)
    open(out_fa, "w").close()
    raise SystemExit(0)

df = pd.read_csv(summary_path, sep="\t")

for col in ["consensus","best_species","pident","qcov","status"]:
    if col not in df.columns:
        raise ValueError(f"Missing column '{col}' in {summary_path}. Found: {list(df.columns)}")

seqs = read_fasta(cons_fa)
cluster_sizes = parse_uc_cluster_sizes(uc_path)

def get_read_count(cons_id):
    if cons_id in cluster_sizes:
        return int(cluster_sizes[cons_id])
    return 0

df["read_count"] = df["consensus"].map(get_read_count).fillna(0).astype(int)
total_reads = int(df["read_count"].sum())
if total_reads > 0:
    df["rel_abundance"] = df["read_count"] / total_reads
else:
    df["rel_abundance"] = 0.0

df["is_rare"] = (df["read_count"] > 0) & (
    (df["read_count"] < min_reads) | (df["rel_abundance"] < min_rel)
)

def is_unclassified(x):
    if pd.isna(x):
        return True
    s = str(x).strip()
    if s == "":
        return True
    low = s.lower()
    return low in {"unclassified", "unknown", "na", "none", "unassigned"}

def novel_row(r):
    pid = float(r["pident"]) if not pd.isna(r["pident"]) else -1.0
    qcov = float(r["qcov"]) if not pd.isna(r["qcov"]) else -1.0
    status = "" if pd.isna(r["status"]) else str(r["status"]).strip().lower()
    unclass = is_unclassified(r["best_species"])
    low_conf = (pid >= 0 and pid < novel_pid) or (qcov >= 0 and qcov < novel_qcov)
    bad_status = status not in {"ok", "pass", "passed", "good"} and status != ""
    if unclass_is_novel and unclass:
        return True
    if low_conf:
        return True
    if bad_status:
        return True
    return False

df["is_novel"] = df.apply(novel_row, axis=1)

reasons = []
for _, r in df.iterrows():
    rs = []
    if bool(r["is_rare"]):
        rs.append("rare")
    if bool(r["is_novel"]):
        status = "" if pd.isna(r["status"]) else str(r["status"]).strip()
        pid = r["pident"]
        qcov = r["qcov"]
        if unclass_is_novel and is_unclassified(r["best_species"]):
            rs.append("unclassified")
        if (not pd.isna(pid)) and float(pid) < novel_pid:
            rs.append(f"low_pid<{novel_pid}")
        if (not pd.isna(qcov)) and float(qcov) < novel_qcov:
            rs.append(f"low_qcov<{novel_qcov}")
        if status and status.strip().lower() not in {"ok","pass","passed","good"}:
            rs.append(f"status={status}")
    reasons.append(";".join(rs))

df["flag_reason"] = reasons

df["blast_url"] = df["consensus"].map(lambda cid: blast_url(seqs.get(cid, "")) if cid in seqs else "")

df.to_csv(out_tsv, sep="\t", index=False)

with open(out_fa, "w") as out:
    for _, r in df.iterrows():
        if not (bool(r["is_rare"]) or bool(r["is_novel"])):
            continue
        cid = str(r["consensus"])
        seq = seqs.get(cid, "")
        if not seq:
            continue
        hdr = f">{cid} reason={r['flag_reason']} best_species={r['best_species']}"
        out.write(hdr + "\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")