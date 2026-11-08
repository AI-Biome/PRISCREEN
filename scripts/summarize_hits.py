#!/usr/bin/env python3
"""Assign a taxon to every consensus sequence from its MMseqs2 hits.

Input: the consensus FASTA of one amplicon and the MMseqs2 ``easy-search`` output
of its database (shared by all amplicons of that database) with the columns
``query,target,theader,pident,qcov,tcov,alnlen,bits,evalue`` (``pident`` is a
percentage, ``qcov`` a fraction). Hits below ``min_pident`` / ``min_qcov`` are
discarded. Of the remaining hits of a query, all within ``top_delta`` of the
best bit score (``bits >= best * (1 - top_delta)``) are the tie band; the taxon
is the lowest rank shared by every hit in the band (species when all agree on
the species). ``status`` is ``PASS_<RANK>``, ``AMBIGUOUS`` (band disagrees
even at the domain), ``BELOW_THRESHOLD`` (hits exist but none passes the
thresholds) or ``NO_HIT``.
"""

import os
import sys

import pandas as pd

COLUMNS = ["query", "target", "theader", "pident", "qcov", "tcov", "alnlen", "bits", "evalue"]
OUT_COLUMNS = ["consensus", "best_species", "assigned_rank", "pident", "qcov", "bits",
               "n_hits_considered", "n_distinct_species", "closest_target", "status"]


def read_fasta_ids(path):
    ids = []
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.append(line[1:].split()[0])
    return ids


def summarize(hits: pd.DataFrame, queries, min_pident: float, min_qcov: float, top_delta: float,
              parse_taxonomy, lowest_common_rank) -> pd.DataFrame:
    """One row per query id (consensus); queries without hits get status NO_HIT."""
    rows = []
    groups = {q: g for q, g in hits.groupby("query", sort=False)} if not hits.empty else {}
    for q in queries:
        dfq = groups.get(q)
        if dfq is not None:
            dfq = dfq.dropna(subset=["pident", "qcov", "bits"]).sort_values(
                ["bits", "pident", "qcov"], ascending=[False, False, False])
        if dfq is None or dfq.empty:
            rows.append(dict(consensus=q, best_species="", assigned_rank="", pident="", qcov="", bits="",
                             n_hits_considered=0, n_distinct_species=0, closest_target="", status="NO_HIT"))
            continue
        closest = dfq.iloc[0]
        passing = dfq[(dfq["pident"] >= min_pident) & (dfq["qcov"] >= min_qcov)]
        if passing.empty:
            rows.append(dict(consensus=q, best_species="", assigned_rank="",
                             pident=float(closest["pident"]), qcov=float(closest["qcov"]), bits=float(closest["bits"]),
                             n_hits_considered=0, n_distinct_species=0, closest_target=str(closest["target"]),
                             status="BELOW_THRESHOLD"))
            continue
        best_bits = float(passing["bits"].max())
        band = passing[passing["bits"] >= best_bits * (1.0 - top_delta)]
        taxa = [parse_taxonomy(h) for h in band["theader"].astype(str)]
        species = sorted({t["species"] for t in taxa if t.get("species")})
        rank, name = lowest_common_rank(taxa)
        top = band.iloc[0]
        rows.append(dict(consensus=q, best_species=name or "", assigned_rank=rank or "",
                         pident=float(top["pident"]), qcov=float(top["qcov"]), bits=float(top["bits"]),
                         n_hits_considered=int(len(band)), n_distinct_species=len(species),
                         closest_target=str(top["target"]),
                         status=f"PASS_{rank.upper()}" if rank else "AMBIGUOUS"))
    return pd.DataFrame(rows, columns=OUT_COLUMNS)


def read_hits(path: str) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=COLUMNS)
    hits = pd.read_csv(path, sep="\t", header=None, names=COLUMNS, quoting=3, dtype=str)
    for c in ["pident", "qcov", "tcov", "alnlen", "bits", "evalue"]:
        hits[c] = pd.to_numeric(hits[c], errors="coerce")
    return hits


def main(smk) -> None:
    sys.path.insert(0, smk.scriptdir)
    from taxonomy import lowest_common_rank, parse_header_taxonomy

    queries = read_fasta_ids(smk.input.cons)
    hits = read_hits(smk.input.hits)
    hits = hits[hits["query"].isin(set(queries))]
    out = summarize(hits, queries, float(smk.params.min_pident), float(smk.params.min_qcov), float(smk.params.top_delta),
                    parse_header_taxonomy, lowest_common_rank)
    out.to_csv(smk.output.tsv, sep="\t", index=False)


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
