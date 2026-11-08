#!/usr/bin/env python3
"""Build the PRISCREEN reference panel from PRISTINE output.

PRISTINE writes the predicted amplicon of every input genome to
``<output>/<group>/primers/amplicons/<locus>_amplicons.fasta`` with headers like
``>Yersinia_pestis_1 accn|CP130901:171501-172828 strand=+ size=1327``.
PRISCREEN needs one record per amplicon variant with the header
``>locus|species=Genus_species|genome=Yersinia_pestis_1``.

The species is the PRISTINE genome id without its trailing ``_<i>`` (the same
convention PRISTINE uses to group genomes into species). Identical sequences of
the same locus and species are written once.

Usage (standalone):
    python scripts/pristine_to_panel.py PRISTINE_OUTPUT_DIR [MORE_DIRS...] -o panel.fasta

The directory may be a PRISTINE output root (several groups), one group
directory, a ``primers/amplicons`` directory, or a lab package with an
``amplicons/`` folder next to ``selected_primers.csv``. The same function is used by
the ``build_panel_from_pristine`` rule of the workflow.

A PRISTINE run that already wrote a ready-made panel (``priscreen_panel.fasta``
at the output root or in a group's ``primers/`` directory, headers
``>locus|species=...|genome=...``) is used as is: that file replaces the
amplicon files of its directory.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

GENOME_SUFFIX = re.compile(r"_\d+$")


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header without '>', sequence) for every record of a FASTA file."""
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


def species_from_genome_id(genome_id: str) -> str:
    """'Yersinia_pestis_12' -> 'Yersinia_pestis'."""
    return GENOME_SUFFIX.sub("", genome_id)


PREBUILT_PANEL = "priscreen_panel.fasta"


def is_prebuilt_panel(path: Path) -> bool:
    return Path(path).name == PREBUILT_PANEL


def find_amplicon_files(roots: Iterable[Path]) -> List[Path]:
    """All <locus>_amplicons.fasta files below the given PRISTINE directories.

    A ``priscreen_panel.fasta`` written by PRISTINE at a root or in a group's ``primers/``
    directory takes the place of the amplicon files it covers."""
    files: List[Path] = []
    for root in roots:
        root = Path(root)
        if root.is_file():
            if not (is_prebuilt_panel(root) or root.name.endswith("_amplicons.fasta")):
                raise ValueError(f"{root}: expected a PRISTINE <locus>_amplicons.fasta or {PREBUILT_PANEL}")
            files.append(root)
            continue
        if not root.is_dir():
            raise FileNotFoundError(f"PRISTINE directory not found: {root}")
        if (root / PREBUILT_PANEL).is_file():          # output root (all groups) or primers/ directory
            files.append(root / PREBUILT_PANEL)
        elif (root / "primers" / PREBUILT_PANEL).is_file():   # one group directory
            files.append(root / "primers" / PREBUILT_PANEL)
        elif root.name == "amplicons":
            files.extend(sorted(root.glob("*_amplicons.fasta")))
        else:
            files.extend(sorted(root.glob("primers/amplicons/*_amplicons.fasta")))
            files.extend(sorted(root.glob("amplicons/*_amplicons.fasta")))       # lab package layout
            for group in sorted(p for p in root.iterdir() if p.is_dir()):
                if (group / "primers" / PREBUILT_PANEL).is_file():
                    files.append(group / "primers" / PREBUILT_PANEL)
                else:
                    files.extend(sorted(group.glob("primers/amplicons/*_amplicons.fasta")))
    uniq = sorted(set(f.resolve() for f in files))
    if not uniq:
        raise FileNotFoundError(
            "no primers/amplicons/*_amplicons.fasta found below " + ", ".join(str(r) for r in roots))
    return uniq


def prebuilt_records(path: Path) -> Iterator[Tuple[str, str, str]]:
    """(locus, species, header) of every record of a ready-made panel; headers are validated."""
    for header, seq in read_fasta(path):
        head = header.split()[0]
        fields = head.split("|")
        species = next((f[len("species="):] for f in fields[1:] if f.startswith("species=")), None)
        if not species or not fields[0]:
            raise ValueError(f"{path}: panel header needs 'locus|species=Genus_species': {header}")
        yield fields[0], species, head, seq


def locus_from_path(path: Path) -> str:
    return path.name[: -len("_amplicons.fasta")]


def panel_records(files: Iterable[Path]) -> List[Tuple[str, str]]:
    """PRISCREEN panel records (header, sequence) from PRISTINE amplicon files."""
    records: List[Tuple[str, str]] = []
    seen: Dict[Tuple[str, str, str], str] = {}
    loci_seen: Dict[str, Path] = {}
    for path in files:
        if is_prebuilt_panel(path):
            for locus, species, head, seq in prebuilt_records(path):
                if locus in loci_seen and loci_seen[locus].parent != path.parent:
                    locus = f"{path.parents[1].name}_{locus}"
                    head = "|".join([locus] + head.split("|")[1:])
                loci_seen[locus] = path
                key = (locus, species, seq.upper())
                if key not in seen:
                    seen[key] = head
                    records.append((head, seq.upper()))
            continue
        locus = locus_from_path(path)
        if "|" in locus:
            raise ValueError(f"locus name must not contain '|': {locus}")
        if locus in loci_seen and loci_seen[locus].parent != path.parent:
            # Same locus name from two PRISTINE groups: prefix with the group directory.
            locus = f"{path.parents[2].name}_{locus}"
        loci_seen[locus] = path
        for header, seq in read_fasta(path):
            genome_id = header.split()[0]
            species = species_from_genome_id(genome_id)
            key = (locus, species, seq.upper())
            if key in seen:
                continue
            seen[key] = genome_id
            records.append((f"{locus}|species={species}|genome={genome_id}", seq.upper()))
    return records


def write_panel(records: Iterable[Tuple[str, str]], out, width: int = 80) -> int:
    n = 0
    for header, seq in records:
        out.write(f">{header}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i:i + width] + "\n")
        n += 1
    return n


def amplicon_names(files: Iterable[Path]) -> List[str]:
    """Amplicon (locus) ids the panel will contain, without reading sequences."""
    return sorted({h.split("|")[0] for h, _ in panel_records(files)})


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pristine_dirs", nargs="+", type=Path, help="PRISTINE output root, group or primers/amplicons directory")
    ap.add_argument("-o", "--output", type=Path, default=None, help="panel FASTA (default: stdout)")
    args = ap.parse_args(argv)
    files = find_amplicon_files(args.pristine_dirs)
    records = panel_records(files)
    if args.output:
        with open(args.output, "w") as out:
            n = write_panel(records, out)
    else:
        n = write_panel(records, sys.stdout)
    loci = {h.split("|")[0] for h, _ in records}
    species = {h.split("|")[1].split("=", 1)[1] for h, _ in records}
    print(f"{n} panel sequences for {len(loci)} amplicon(s) and {len(species)} species from {len(files)} PRISTINE file(s)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
