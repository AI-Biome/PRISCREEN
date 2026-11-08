#!/usr/bin/env python3
"""Parse taxonomy from reference database headers.

PRISCREEN identifies consensus sequences against two kinds of databases and
needs a lineage for every hit:

* the PRISTINE panel: ``locus|species=Genus_species|genome=...``
* rRNA databases with the lineage in the header description, e.g.
  SILVA ``>AB001.1.1500 Bacteria;Pseudomonadota;Gammaproteobacteria;...;Escherichia coli``,
  GTDB SSU ``>RS_GCF_000... d__Bacteria;p__...;s__Escherichia coli [location=...]``,
  SINTAX-style ``>id;tax=d:Bacteria,p:...,s:Escherichia_coli;`` and
  ``>id|tax=Bacteria;...;Escherichia coli`` or ``>id|genus=X|species=Y``.

A lineage is a list of (rank, name) from the broadest to the finest level.
Ranks are the seven canonical ones; levels that cannot be named (SILVA has a
variable number of levels) are called ``level<i>``.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Sequence, Tuple

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_ALIASES = {"kingdom": "domain", "superkingdom": "domain", "d": "domain", "k": "domain",
                "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
UNCLASSIFIED = re.compile(r"^(unclassified|uncultured|unknown|unidentified|unassigned|na|none|metagenome|environmental)\b", re.I)

Lineage = List[Tuple[str, str]]


def rank_index(rank: Optional[str]) -> int:
    """Position of a rank in RANKS, or -1 for unknown ranks."""
    if rank is None:
        return -1
    r = RANK_ALIASES.get(str(rank).strip().lower(), str(rank).strip().lower())
    return RANKS.index(r) if r in RANKS else -1


def _clean(name: str) -> Optional[str]:
    name = name.strip().strip('"').strip().rstrip(";").strip()
    if not name or UNCLASSIFIED.match(name):
        return None
    return name


def _prefixed(parts: Sequence[str]) -> Optional[Lineage]:
    """Lineage from 'd__Bacteria' / 'd:Bacteria' style parts; None if not prefixed."""
    out: Lineage = []
    for p in parts:
        m = re.match(r"^\s*([a-zA-Z]{1,2})(?:__|:)(.*)$", p)
        if not m:
            return None
        rank = RANK_ALIASES.get(m.group(1).lower(), m.group(1).lower())
        name = _clean(m.group(2).replace("_", " "))
        out.append((rank, name) if name else (rank, None))
    return out if out else None


def _positional(parts: Sequence[str]) -> Lineage:
    """Unprefixed lineage (SILVA has no rank labels).

    Seven levels are read as domain..species. Longer lineages (SILVA eukaryotes,
    bacterial subclasses) keep the first level as domain, the last two as genus
    and species, and the levels in between unnamed (``level<i>``), so that the
    lowest common rank of such hits is still found by position.
    """
    names = [_clean(p) for p in parts]
    out: Lineage = []
    if len(names) == 7:
        for i, n in enumerate(names):
            out.append((RANKS[i], n))
    elif len(names) == 1:
        out.append(("species", names[0]))
    elif len(names) < 7:
        # SILVA: "<path>;<organism name>" with a variable path depth. The last element is the
        # organism (species), the first the domain; the rest fill the ranks from the top.
        out.append(("domain", names[0]))
        for i, n in enumerate(names[1:-1], start=1):
            out.append((RANKS[i], n))
        out.append(("species", names[-1]))
    else:
        out.append(("domain", names[0]))
        for i, n in enumerate(names[1:-2], start=2):
            out.append((f"level{i}", n))
        out.append(("genus", names[-2]))
        out.append(("species", names[-1]))
    return out


def parse_header_taxonomy(header: str) -> Dict[str, Optional[str]]:
    """Map rank -> name for one database header (without '>').

    Unknown names are None. Always contains the keys in RANKS; ``species`` is
    the finest resolved level, ``genus`` is inferred from a binomial species name.
    """
    s = str(header).strip()
    tax: Dict[str, Optional[str]] = {r: None for r in RANKS}
    lineage: Optional[Lineage] = None

    ident = s.split()[0] if s else ""
    m = re.search(r"(?:^|[|;\s])tax=([^|;]+(?:;[^|;]*)*)", s)
    if m and ":" in m.group(1) and "," in m.group(1):          # SINTAX: tax=d:Bacteria,p:...
        lineage = _prefixed([p for p in m.group(1).split(",") if p.strip()])
    elif "|" in ident and "=" in ident.split("|", 1)[1]:        # PRISCREEN panel: locus|species=X|... (the id token)
        kv = {}
        for field in ident.split("|")[1:]:
            if "=" in field:
                k, v = field.split("=", 1)
                kv[k.strip().lower()] = v.strip()
        if "tax" in kv:
            parts = [p for p in kv["tax"].split(";") if p.strip()]
            lineage = _prefixed(parts) or _positional(parts)
        else:
            for k, v in kv.items():
                r = RANK_ALIASES.get(k, k)
                if r in RANKS:
                    tax[r] = _clean(v.replace("_", " "))
    else:                                                        # SILVA / GTDB: "id Lineage;..." or "id;Lineage"
        desc = s.split(None, 1)[1] if " " in s else (s.split(";", 1)[1] if ";" in s else "")
        desc = re.sub(r"\s*\[.*?\]\s*", " ", desc).strip()     # GTDB [location=...] annotations
        parts = [p for p in desc.split(";") if p.strip()]
        if parts:
            lineage = _prefixed(parts) or _positional(parts)

    if lineage:
        for rank, name in lineage:
            if rank in tax and name:
                tax[rank] = name

    if tax["species"] and tax["genus"] is None:
        tokens = re.split(r"[\s_]+", tax["species"].strip())
        if len(tokens) >= 2 and tokens[0][0].isupper():        # binomial: Genus species
            tax["genus"] = tokens[0]
    if all(v is None for v in tax.values()):
        tax["species"] = s.split()[0] if s else None                     # bare id, e.g. a plain FASTA
    return tax


def lowest_common_rank(taxa: Sequence[Dict[str, Optional[str]]]) -> Tuple[Optional[str], Optional[str]]:
    """(rank, name) of the finest rank on which all taxa agree, or (None, None).

    A rank that is unknown for every taxon (the panel only names species) is
    skipped; a rank that is known for some taxa and unknown for others, or on
    which the taxa disagree, ends the search.
    """
    deepest: Tuple[Optional[str], Optional[str]] = (None, None)
    for rank in RANKS:
        values = {t.get(rank) for t in taxa}
        if values == {None}:
            continue
        if len(values) != 1 or None in values:
            break
        deepest = (rank, values.pop())
    return deepest
