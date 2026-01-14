#!/usr/bin/env python3
import os
import re
import pandas as pd

hits_path = snakemake.input.hits
out_tsv = snakemake.output.tsv

min_pid = float(snakemake.params.min_pid)
min_qcov = float(snakemake.params.min_qcov)
top_delta = float(snakemake.params.top_delta)

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

def parse_target_taxonomy(target: str):
    tax = {r: None for r in RANKS}
    if target is None:
        return tax

    s = str(target).strip()

    if "tax=" in s:
        m = re.search(r"(?:^|\|)tax=([^|]+)", s)
        if m:
            parts = [p.strip() for p in m.group(1).split(";") if p.strip()]
            for i, p in enumerate(parts):
                if i >= len(RANKS):
                    break
                tax[RANKS[i]] = p
            return tax

    fields = s.split("|")
    kv = {}
    for f in fields[1:]:
        if "=" in f:
            k, v = f.split("=", 1)
            kv[k.strip().lower()] = v.strip()

    for r in RANKS:
        if r in kv:
            tax[r] = kv[r]

    if tax["species"] is None and "species" in kv:
        tax["species"] = kv["species"]

    if tax["species"] is None:
        tax["species"] = s

    if tax["genus"] is None and tax["species"]:
        sp = str(tax["species"]).strip()
        if sp:
            g = re.split(r"[_\s]+", sp)[0]
            tax["genus"] = g if g else None

    return tax

def lca_rank_and_name(tax_list):
    lineages = []
    for t in tax_list:
        lineage = []
        for r in RANKS:
            v = t.get(r)
            lineage.append(v if v not in (None, "", "NA", "na", "None") else None)
        lineages.append(lineage)

    deepest = None
    for i, r in enumerate(RANKS):
        vals = [lin[i] for lin in lineages]
        if any(v is None for v in vals):
            break
        if len(set(vals)) == 1:
            deepest = (r, vals[0])
        else:
            break

    return deepest

def status_from_rank(rank):
    if rank is None:
        return "AMBIGUOUS"
    return f"PASS_{rank.upper()}"

if (not os.path.exists(hits_path)) or os.path.getsize(hits_path) == 0:
    pd.DataFrame(columns=[
        "consensus", "best_species", "assigned_rank",
        "pident", "qcov", "bits",
        "n_hits_considered", "n_distinct_species",
        "status"
    ]).to_csv(out_tsv, sep="\t", index=False)
    raise SystemExit(0)

hits = pd.read_csv(
    hits_path, sep="\t", header=None,
    names=["query", "target", "pident", "qcov", "bits"]
)

for c in ["pident", "qcov", "bits"]:
    hits[c] = pd.to_numeric(hits[c], errors="coerce")

rows = []
for q, dfq in hits.groupby("query", sort=False):
    dfq = dfq.dropna(subset=["pident", "qcov", "bits"])
    dfq = dfq[(dfq["pident"] >= min_pid) & (dfq["qcov"] >= min_qcov)]

    if dfq.empty:
        rows.append({
            "consensus": q,
            "best_species": "",
            "assigned_rank": "",
            "pident": "",
            "qcov": "",
            "bits": "",
            "n_hits_considered": 0,
            "n_distinct_species": 0,
            "status": "NO_HIT"
        })
        continue

    best_bits = float(dfq["bits"].max())
    band = dfq[dfq["bits"] >= (best_bits - top_delta)].copy()
    band = band.sort_values(["bits", "pident", "qcov"], ascending=[False, False, False])

    tax_list = [parse_target_taxonomy(t) for t in band["target"].tolist()]
    species_list = [t.get("species") for t in tax_list if t.get("species") not in (None, "", "NA", "na", "None")]
    distinct_species = sorted(set(species_list))

    top = band.iloc[0]
    assigned_rank = None
    assigned_name = None

    if len(distinct_species) == 1:
        assigned_rank = "species"
        assigned_name = distinct_species[0]
    else:
        lca = lca_rank_and_name(tax_list)
        if lca is not None:
            assigned_rank, assigned_name = lca
        else:
            assigned_rank, assigned_name = None, None

    rows.append({
        "consensus": q,
        "best_species": assigned_name if assigned_name is not None else "",
        "assigned_rank": assigned_rank if assigned_rank is not None else "",
        "pident": float(top["pident"]),
        "qcov": float(top["qcov"]),
        "bits": float(top["bits"]),
        "n_hits_considered": int(band.shape[0]),
        "n_distinct_species": int(len(distinct_species)),
        "status": status_from_rank(assigned_rank)
    })

out = pd.DataFrame(rows, columns=[
    "consensus", "best_species", "assigned_rank",
    "pident", "qcov", "bits",
    "n_hits_considered", "n_distinct_species",
    "status"
])
out.to_csv(out_tsv, sep="\t", index=False)
