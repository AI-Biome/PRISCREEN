#!/usr/bin/env python3

import os
import urllib.parse
import pandas as pd

# Per-amplicon summary table produced by summarize_hits.py
summary_path = snakemake.input.summary

# Consensus FASTA for this sample+amplicon
cons_fa = snakemake.input.cons

# VSEARCH .uc clustering file, used to estimate read support per consensus
uc_path = snakemake.input.uc

# Output files
out_tsv = snakemake.output.tsv
out_fa = snakemake.output.fasta

# Minimum absolute number of reads required not to be considered rare
min_reads = int(snakemake.params.min_reads)

# Minimum relative abundance required not to be considered rare
min_rel = float(snakemake.params.min_rel)

# Identity threshold below which a sequence may be considered potentially novel
novel_pid = float(snakemake.params.novel_pid)

# Query coverage threshold below which a sequence may be considered potentially novel
novel_qcov = float(snakemake.params.novel_qcov)

# Whether "unclassified" assignments should count as potentially novel
treat_unclassified_as_novel = bool(snakemake.params.treat_unclassified_as_novel)

# Lowest acceptable resolved rank before a result is flagged as low resolution
# Example: if "genus", then family/order/class... will be flagged as low resolution
min_resolved_rank = str(snakemake.params.min_resolved_rank)

# Taxonomic rank ordering used for low-resolution detection

RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def rank_index(r):
    """
    Convert a rank string to its numeric position in RANK_ORDER.

    Returns:
        -1 if the rank is missing or invalid
        0..6 for valid ranks
    """
    if r is None:
        return -1
    s = str(r).strip().lower()
    if s not in RANK_ORDER:
        return -1
    return RANK_ORDER.index(s)


# Validate the configured minimum resolved rank early
if rank_index(min_resolved_rank) < 0:
    raise ValueError(
        f"Invalid min_resolved_rank='{min_resolved_rank}'. Allowed: {RANK_ORDER}"
    )


def read_fasta(path):
    """
    Read a FASTA file into a dictionary: {sequence_id: sequence}.
    """
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
                # Save the previous record before starting a new one
                if name is not None:
                    seqs[name] = "".join(parts)

                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

        # Save the final record
        if name is not None:
            seqs[name] = "".join(parts)

    return seqs


def parse_uc_cluster_sizes(path):
    """
    Parse a VSEARCH .uc file and estimate cluster sizes.
    """
    sizes = {}

    if not path or (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return sizes

    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")

            # Only "H" lines represent hits assigned to a centroid
            if cols[0] != "H":
                continue

            # In .uc format, column 10 (0-based index 9) is the centroid label
            centroid = cols[9]
            sizes[centroid] = sizes.get(centroid, 0) + 1

    # Add 1 for the centroid sequence itself
    for c in list(sizes.keys()):
        sizes[c] = sizes[c] + 1

    return sizes


def blast_url(seq):
    """
    Build a direct NCBI BLAST URL for a given nucleotide sequence.
    """
    q = urllib.parse.quote(seq)
    return (
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"
        "PAGE_TYPE=BlastSearch&PROGRAM=blastn&QUERY=" + q
    )


def is_unclassified(x):
    """
    Return True if a taxonomic label should be treated as unclassified/unknown.
    """
    if pd.isna(x):
        return True

    s = str(x).strip()
    if s == "":
        return True

    low = s.lower()
    return low in {"unclassified", "unknown", "na", "none", "unassigned"}


# Expected base columns from summarize_hits.py
base_cols = [
    "consensus",
    "best_species",
    "assigned_rank",
    "pident",
    "qcov",
    "bits",
    "status"
]

# ---------------------------------------------------------------------
# Handle empty input summary
# ---------------------------------------------------------------------

# If the summary file is missing or empty, create empty outputs and exit cleanly
if (not os.path.exists(summary_path)) or os.path.getsize(summary_path) == 0:
    pd.DataFrame(columns=base_cols + [
        "read_count",
        "rel_abundance",
        "is_rare",
        "is_low_resolution",
        "is_novel",
        "flag_reason",
        "blast_url"
    ]).to_csv(out_tsv, sep="\t", index=False)

    open(out_fa, "w").close()
    raise SystemExit(0)

# ---------------------------------------------------------------------
# Load summary table
# ---------------------------------------------------------------------

df = pd.read_csv(summary_path, sep="\t")

# Check that required columns are present
missing = [c for c in ["consensus", "best_species", "pident", "qcov", "status"] if c not in df.columns]
if missing:
    raise ValueError(
        f"Missing columns in {summary_path}: {missing}. Found: {list(df.columns)}"
    )

# Add optional columns if they are not present
if "assigned_rank" not in df.columns:
    df["assigned_rank"] = ""
if "bits" not in df.columns:
    df["bits"] = ""
if "status" not in df.columns:
    df["status"] = ""

# ---------------------------------------------------------------------
# Load auxiliary data
# ---------------------------------------------------------------------

# Consensus sequences
seqs = read_fasta(cons_fa)

# Read support derived from clustering
cluster_sizes = parse_uc_cluster_sizes(uc_path)


def get_read_count(cons_id):
    """
    Return the number of reads supporting a given consensus ID.
    """
    if cons_id in cluster_sizes:
        return int(cluster_sizes[cons_id])
    return 0


# ---------------------------------------------------------------------
# Compute read support and relative abundance
# ---------------------------------------------------------------------

# Per-consensus absolute support
df["read_count"] = df["consensus"].map(get_read_count).fillna(0).astype(int)

# Total reads represented by all consensuses in this sample+amplicon
total_reads = int(df["read_count"].sum())

# Relative abundance within this sample+amplicon
if total_reads > 0:
    df["rel_abundance"] = df["read_count"] / total_reads
else:
    df["rel_abundance"] = 0.0

# ---------------------------------------------------------------------
# Rare taxa flag
# ---------------------------------------------------------------------

# A taxon is considered rare if:
#   - it has at least one supporting read
#   - and either low absolute support or low relative abundance
df["is_rare"] = (df["read_count"] > 0) & (
    (df["read_count"] < min_reads) | (df["rel_abundance"] < min_rel)
)

# ---------------------------------------------------------------------
# Low-resolution flag
# ---------------------------------------------------------------------

# A taxon is low resolution if it is resolved only above the configured threshold.
# Example with min_resolved_rank="genus":
#   family/order/class/phylum/kingdom -> low resolution
#   genus/species -> not low resolution
df["is_low_resolution"] = df["assigned_rank"].apply(
    lambda r: rank_index(r) >= 0 and rank_index(r) < rank_index(min_resolved_rank)
)

# ---------------------------------------------------------------------
# Novelty flag
# ---------------------------------------------------------------------

def novel_row(r):
    """
    Decide whether a row should be flagged as potentially novel.

    Criteria:
      - NO_HIT or AMBIGUOUS status
      - unclassified taxonomy (if configured to count as novel)
      - low pident or low qcov
      - missing pident/qcov
    """
    status = "" if pd.isna(r["status"]) else str(r["status"]).strip().upper()

    if status in {"NO_HIT", "AMBIGUOUS"}:
        return True

    if treat_unclassified_as_novel and is_unclassified(r["best_species"]):
        return True

    try:
        pid = float(r["pident"])
    except Exception:
        pid = None

    try:
        qcov = float(r["qcov"])
    except Exception:
        qcov = None

    if pid is None or qcov is None:
        return True

    if pid < novel_pid:
        return True

    if qcov < novel_qcov:
        return True

    return False


df["is_novel"] = df.apply(novel_row, axis=1)

# ---------------------------------------------------------------------
# Human-readable reasons for flagging
# ---------------------------------------------------------------------

reasons = []

for _, r in df.iterrows():
    rs = []

    # Rare
    if bool(r["is_rare"]):
        rs.append("rare")

    # Low taxonomic resolution
    if bool(r["is_low_resolution"]):
        rs.append(f"low_resolution<{min_resolved_rank}")

    # Potential novelty / uncertain identification
    if bool(r["is_novel"]):
        status = "" if pd.isna(r["status"]) else str(r["status"]).strip().upper()

        if status in {"NO_HIT", "AMBIGUOUS"}:
            rs.append(status.lower())

        if treat_unclassified_as_novel and is_unclassified(r["best_species"]):
            rs.append("unclassified")

        try:
            pid = float(r["pident"])
            if pid < novel_pid:
                rs.append(f"low_pid<{novel_pid}")
        except Exception:
            rs.append("pid_missing")

        try:
            qcov = float(r["qcov"])
            if qcov < novel_qcov:
                rs.append(f"low_qcov<{novel_qcov}")
        except Exception:
            rs.append("qcov_missing")

    reasons.append(";".join(rs))

df["flag_reason"] = reasons

# ---------------------------------------------------------------------
# BLAST links
# ---------------------------------------------------------------------

# Build a direct BLAST URL for each consensus sequence
df["blast_url"] = df["consensus"].map(
    lambda cid: blast_url(seqs.get(cid, "")) if cid in seqs and seqs.get(cid, "") else ""
)

# ---------------------------------------------------------------------
# Write tabular output
# ---------------------------------------------------------------------

df.to_csv(out_tsv, sep="\t", index=False)

# ---------------------------------------------------------------------
# Write FASTA of flagged sequences only
# ---------------------------------------------------------------------

with open(out_fa, "w") as out:
    for _, r in df.iterrows():
        # Include any sequence flagged as rare, low resolution, or novel
        if not (bool(r["is_rare"]) or bool(r["is_low_resolution"]) or bool(r["is_novel"])):
            continue

        cid = str(r["consensus"])
        seq = seqs.get(cid, "")

        if not seq:
            continue

        hdr = f">{cid} reason={r['flag_reason']} best_species={r['best_species']}"
        out.write(hdr + "\n")

        # Wrap sequence to 80 characters per line
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")
