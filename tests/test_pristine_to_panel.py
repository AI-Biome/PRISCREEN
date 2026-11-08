import subprocess
import sys

import pytest
from pathlib import Path

from pristine_to_panel import amplicon_names, find_amplicon_files, panel_records, species_from_genome_id


def test_species_from_genome_id():
    assert species_from_genome_id("Yersinia_pestis_12") == "Yersinia_pestis"
    assert species_from_genome_id("Yersinia_pestis") == "Yersinia_pestis"


def test_find_files_from_root_group_and_amplicon_dir(pristine_toy):
    files = find_amplicon_files([pristine_toy])
    assert len(files) == 3
    assert len(find_amplicon_files([pristine_toy / "Yersinia_pseudotuberculosis"])) == 2
    assert len(find_amplicon_files([pristine_toy / "Yersinia_pseudotuberculosis" / "primers" / "amplicons"])) == 2


def test_panel_records_dedup_and_headers(pristine_toy):
    recs = panel_records(find_amplicon_files([pristine_toy]))
    headers = [h for h, _ in recs]
    # identical sequences of one species and locus are written once
    assert headers.count("g1_locusA|species=Yersinia_pseudotuberculosis|genome=Yersinia_pseudotuberculosis_1") == 1
    assert not any("Yersinia_pseudotuberculosis_2" in h for h in headers)
    assert "g1_locusA|species=Yersinia_pestis|genome=Yersinia_pestis_1" in headers
    assert "g2_locusC|species=Escherichia_coli|genome=Escherichia_coli_1" in headers
    assert all("|species=" in h for h in headers)
    assert amplicon_names(find_amplicon_files([pristine_toy])) == ["g1_locusA", "g1_locusB", "g2_locusC"]


def test_cli(pristine_toy, tmp_path, root):
    out = tmp_path / "panel.fasta"
    r = subprocess.run([sys.executable, str(root / "scripts" / "pristine_to_panel.py"), str(pristine_toy), "-o", str(out)],
                       capture_output=True, text=True, check=True)
    assert "3 amplicon(s)" in r.stderr
    text = out.read_text()
    headers = [l for l in text.splitlines() if l.startswith(">")]
    assert len(headers) == 6 and all("|species=" in h for h in headers)
    assert any(h.startswith(">g1_locusA|species=Yersinia_pseudotuberculosis|") for h in headers)


def test_prebuilt_panel_replaces_amplicon_files(pristine_toy):
    """A priscreen_panel.fasta written by PRISTINE is used as is, for its group only."""
    from pristine_to_panel import find_amplicon_files, panel_records
    group = pristine_toy / "Yersinia_pseudotuberculosis"
    (group / "primers" / "priscreen_panel.fasta").write_text(
        ">g1_locusA|species=Yersinia_pestis|genome=Yersinia_pestis_1|role=non-target\nACGA\n"
        ">g1_locusA|species=Yersinia_pseudotuberculosis|genome=Yersinia_pseudotuberculosis_1|role=target\nACGT\n"
        ">g1_locusA|species=Yersinia_pseudotuberculosis|genome=Yersinia_pseudotuberculosis_2|role=target\nACGT\n")
    files = find_amplicon_files([pristine_toy])
    assert [f.name for f in files] == ["g2_locusC_amplicons.fasta", "priscreen_panel.fasta"]
    recs = panel_records(files)
    heads = [h for h, _ in recs]
    assert heads == ["g2_locusC|species=Salmonella_enterica|genome=Salmonella_enterica_1",
                     "g2_locusC|species=Escherichia_coli|genome=Escherichia_coli_1",
                     "g1_locusA|species=Yersinia_pestis|genome=Yersinia_pestis_1|role=non-target",
                     "g1_locusA|species=Yersinia_pseudotuberculosis|genome=Yersinia_pseudotuberculosis_1|role=target"]
    # a root-level panel covers everything below it
    (pristine_toy / "priscreen_panel.fasta").write_text(">x|species=Only_one|genome=Only_one_1\nAAAA\n")
    assert [h for h, _ in panel_records(find_amplicon_files([pristine_toy]))] == ["x|species=Only_one|genome=Only_one_1"]
    (pristine_toy / "priscreen_panel.fasta").write_text(">x|genome=Only_one_1\nAAAA\n")
    with pytest.raises(ValueError):
        panel_records(find_amplicon_files([pristine_toy]))


def test_lab_package_layout(tmp_path):
    """A lab package has amplicons/ directly below the root (no primers/)."""
    from pristine_to_panel import find_amplicon_files, panel_records
    amp = tmp_path / "amplicons"
    amp.mkdir()
    (amp / "locusZ_amplicons.fasta").write_text(">Mycobacterium_avium_paratuberculosis_1 NC_1:1-8 strand=+ size=8\nACGTACGT\n")
    recs = panel_records(find_amplicon_files([tmp_path]))
    assert recs == [("locusZ|species=Mycobacterium_avium_paratuberculosis|genome=Mycobacterium_avium_paratuberculosis_1", "ACGTACGT")]
