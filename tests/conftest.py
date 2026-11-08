import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="session")
def root() -> Path:
    return ROOT


@pytest.fixture
def pristine_toy(tmp_path: Path) -> Path:
    """A fake PRISTINE output tree with two groups and PRISTINE-style amplicon headers."""
    seqs = {
        "g1_locusA": [("Yersinia_pseudotuberculosis_1", "ACGT" * 50), ("Yersinia_pseudotuberculosis_2", "ACGT" * 50),
                      ("Yersinia_pestis_1", "ACGA" * 50)],
        "g1_locusB": [("Yersinia_pseudotuberculosis_1", "TTGA" * 60), ("Yersinia_pestis_3", "TTGC" * 60)],
    }
    out = tmp_path / "pristine_out"
    amp = out / "Yersinia_pseudotuberculosis" / "primers" / "amplicons"
    amp.mkdir(parents=True)
    for locus, recs in seqs.items():
        with open(amp / f"{locus}_amplicons.fasta", "w") as fh:
            for gid, s in recs:
                fh.write(f">{gid} accn|CP1:100-{100 + len(s)} strand=+ size={len(s)}\n{s}\n")
    amp2 = out / "Salmonella_enterica" / "primers" / "amplicons"
    amp2.mkdir(parents=True)
    with open(amp2 / "g2_locusC_amplicons.fasta", "w") as fh:
        fh.write(">Salmonella_enterica_1 accn|X:1-200 strand=- size=200\n" + "GATC" * 50 + "\n")
        fh.write(">Escherichia_coli_1 accn|Y:1-200 strand=+ size=200\n" + "GATG" * 50 + "\n")
    return out
