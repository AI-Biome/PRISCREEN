#!/usr/bin/env python3
"""Flag rare, poorly resolved and potentially novel consensus sequences.

Reads the per-amplicon ``summary.tsv`` of summarize_hits.py, the consensus
FASTA and ``cluster_sizes.tsv`` (consensus id -> cluster size and reads mapped)
and adds:

* ``read_count`` / ``cluster_size`` / ``rel_abundance``: reads mapped to the
  consensus in the final polishing alignment, reads in its vsearch cluster, and
  the share of the mapped reads of this amplicon in this sample;
* ``is_rare``: fewer than ``min_reads`` reads or below ``min_rel_abundance``;
* ``is_low_resolution``: assigned above ``min_resolved_rank`` (e.g. only to a family);
* ``is_novel``: no hit, below the identification thresholds, ambiguous, or
  identity/coverage to the closest reference below ``novel_pident`` (percent) /
  ``novel_qcov`` (fraction);
* ``separated_by`` (snp, indel, snp+indel) for haplotypes of a split cluster; a
  haplotype separated by indels only gets the ``flag_reason`` ``indel_only_split``
  (ONT reads are least reliable there) without counting as novel;
* ``flag_reason`` and a ready-made NCBI BLAST link.

Flagged sequences are also written to a FASTA file for manual follow-up.
"""

import os
import sys
import urllib.parse

import pandas as pd

BASE_COLUMNS = ["consensus", "best_species", "assigned_rank", "pident", "qcov", "bits",
                "n_hits_considered", "n_distinct_species", "closest_target", "status"]
FLAG_COLUMNS = ["read_count", "cluster_size", "separated_by", "rel_abundance", "is_rare", "is_low_resolution", "is_novel",
                "flag_reason", "blast_url"]


def read_fasta(path: str) -> dict:
    seqs, name, parts = {}, None, []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return seqs
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def read_cluster_sizes(path: str):
    """{consensus: (cluster_size, reads_mapped, separated_by)} from cluster_sizes.tsv."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    df = pd.read_csv(path, sep="\t", dtype={"consensus": str, "separated_by": str}, keep_default_na=False)
    reads = df["reads_mapped"] if "reads_mapped" in df.columns else df["cluster_size"]
    sep = df["separated_by"] if "separated_by" in df.columns else [""] * len(df)
    return {c: (int(s), int(r), str(x)) for c, s, r, x in zip(df["consensus"], df["cluster_size"], reads, sep)}


def blast_url(seq: str) -> str:
    return "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn&QUERY=" + urllib.parse.quote(seq)


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def flag(df: pd.DataFrame, sizes: dict, min_reads: int, min_rel: float, novel_pid: float, novel_qcov: float,
         min_resolved_rank: str, rank_index) -> pd.DataFrame:
    """Add the flag columns; the rare test is skipped when no read counts are available (sizes empty)."""
    df = df.copy()
    for c in BASE_COLUMNS:
        if c not in df.columns:
            df[c] = ""
    df["read_count"] = df["consensus"].map(lambda c: int(sizes.get(str(c), (0, 0, ""))[1]))
    df["cluster_size"] = df["consensus"].map(lambda c: int(sizes.get(str(c), (0, 0, ""))[0]))
    df["separated_by"] = df["consensus"].map(lambda c: str(sizes.get(str(c), (0, 0, ""))[2]))
    total = int(df["read_count"].sum())
    df["rel_abundance"] = df["read_count"] / total if total > 0 else 0.0
    if sizes:
        df["is_rare"] = (df["read_count"] < min_reads) | (df["rel_abundance"] < min_rel)
    else:
        df["is_rare"] = False
    min_rank_i = rank_index(min_resolved_rank)
    df["is_low_resolution"] = df["assigned_rank"].map(lambda r: 0 <= rank_index(r) < min_rank_i)

    reasons, novel = [], []
    for _, r in df.iterrows():
        status = str(r["status"]).strip().upper() if not pd.isna(r["status"]) else ""
        rs = []
        if bool(r["is_rare"]):
            rs.append("rare")
        if bool(r["is_low_resolution"]):
            rs.append(f"low_resolution<{min_resolved_rank}")
        if str(r["separated_by"]) == "indel":
            rs.append("indel_only_split")
        nv = []
        if status in {"NO_HIT", "AMBIGUOUS", "BELOW_THRESHOLD"}:
            nv.append(status.lower())
        pid, qcov = to_float(r["pident"]), to_float(r["qcov"])
        if status != "NO_HIT":
            if pid is None:
                nv.append("pid_missing")
            elif pid < novel_pid:
                nv.append(f"low_pid<{novel_pid}")
            if qcov is None:
                nv.append("qcov_missing")
            elif qcov < novel_qcov:
                nv.append(f"low_qcov<{novel_qcov}")
        novel.append(bool(nv))
        reasons.append(";".join(rs + nv))
    df["is_novel"] = novel
    df["flag_reason"] = reasons
    return df


def main(smk) -> None:
    sys.path.insert(0, smk.scriptdir)
    from taxonomy import RANKS, rank_index

    min_resolved_rank = str(smk.params.min_resolved_rank)
    if rank_index(min_resolved_rank) < 0:
        raise ValueError(f"novelty.min_resolved_rank={min_resolved_rank!r} is not one of {RANKS}")

    summary_path, out_tsv, out_fa = smk.input.summary, smk.output.tsv, smk.output.fasta
    if not os.path.exists(summary_path) or os.path.getsize(summary_path) == 0:
        pd.DataFrame(columns=BASE_COLUMNS + FLAG_COLUMNS).to_csv(out_tsv, sep="\t", index=False)
        open(out_fa, "w").close()
        return
    df = pd.read_csv(summary_path, sep="\t")
    if df.empty:
        pd.DataFrame(columns=BASE_COLUMNS + FLAG_COLUMNS).to_csv(out_tsv, sep="\t", index=False)
        open(out_fa, "w").close()
        return
    seqs = read_fasta(smk.input.cons)
    sizes = read_cluster_sizes(smk.input.sizes)
    df = flag(df, sizes, int(smk.params.min_reads), float(smk.params.min_rel), float(smk.params.novel_pid),
              float(smk.params.novel_qcov), min_resolved_rank, rank_index)
    df["blast_url"] = df["consensus"].map(lambda c: blast_url(seqs[c]) if seqs.get(c) else "")
    df.to_csv(out_tsv, sep="\t", index=False)
    with open(out_fa, "w") as out:
        for _, r in df.iterrows():
            if not (r["is_rare"] or r["is_low_resolution"] or r["is_novel"]):
                continue
            seq = seqs.get(str(r["consensus"]), "")
            if not seq:
                continue
            out.write(f">{r['consensus']} reason={r['flag_reason']} best_species={r['best_species']} reads={r['read_count']}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
