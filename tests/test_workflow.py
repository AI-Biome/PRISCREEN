"""Workflow tests: parsing/dry run always; the full toy run only with PRISCREEN_FULL_TEST=1 (needs ~5 min)."""
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest


def snakemake(root, workdir, *args):
    cmd = ["snakemake", "--snakefile", str(root / "Snakefile"), "--directory", str(workdir), *args]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=root)


def test_dry_run(root, tmp_path):
    r = snakemake(root, tmp_path, "-n", "--cores", "2")
    assert r.returncode == 0, r.stderr + r.stdout
    assert "medaka" in r.stdout and "racon" in r.stdout and "read_votes" in r.stdout and "mmseqs_search" in r.stdout
    assert "haplotype_split" in r.stdout


def test_dry_run_with_container(root, tmp_path):
    r = snakemake(root, tmp_path, "-n", "--cores", "2", "--profile", str(root / "profiles" / "local"))
    assert r.returncode == 0, r.stderr + r.stdout
    assert "quay.io/biocontainers/medaka" in r.stdout + r.stderr and "quay.io/biocontainers/racon" in r.stdout + r.stderr


def test_dry_run_without_haplotypes(root, tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text("haplotypes:\n  enable: false\n")
    r = snakemake(root, tmp_path, "-n", "--cores", "2", "--configfile", str(cfg))
    assert r.returncode == 0, r.stderr + r.stdout
    assert "haplotype_split" not in r.stdout and "summarize_hits" in r.stdout


def test_dry_run_from_pristine_output(root, tmp_path, pristine_toy):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(f"""
samples: {root / 'config' / 'samples.tsv'}
panel:
  pristine_output: {pristine_toy}
identify:
  enable: true
""")
    r = snakemake(root, tmp_path, "-n", "--cores", "2", "--configfile", str(cfg))
    assert r.returncode == 0, r.stderr + r.stdout
    assert "build_panel_from_pristine" in r.stdout
    assert "g1_locusA" in r.stdout and "g2_locusC" in r.stdout


@pytest.mark.skipif(os.environ.get("PRISCREEN_FULL_TEST") != "1", reason="set PRISCREEN_FULL_TEST=1 for the full toy run")
def test_full_toy_run(root, tmp_path):
    r = snakemake(root, tmp_path, "--cores", "4")
    assert r.returncode == 0, r.stderr + r.stdout
    summary = pd.read_csv(tmp_path / "results/identify/Test2/species_summary.tsv", sep="\t")
    amp1 = summary.loc[(summary.amplicon == "amp1") & (summary.status == "PASS_SPECIES")]
    assert set(amp1.best_species) == {"SpeciesA", "SpeciesB", "SpeciesD"}
    # SpeciesD differs from SpeciesA by three indels only and shares its cluster: separated by the haplotype step
    sizes1 = pd.read_csv(tmp_path / "results/haplotypes/Test2/amp1/cluster_sizes.tsv", sep="\t", keep_default_na=False).set_index("consensus")
    d_cons = amp1.loc[amp1.best_species == "SpeciesD", "consensus"].iloc[0]
    assert "_h" in d_cons and "indel" in sizes1.loc[d_cons, "separated_by"]
    assert abs(sizes1.loc[d_cons, "reads_mapped"] - 30) <= 3
    novelty2 = pd.read_csv(tmp_path / "results/identify/Test2/novelty_summary.tsv", sep="\t", keep_default_na=False)
    assert "indel_only_split" in novelty2.loc[novelty2.consensus == d_cons, "flag_reason"].iloc[0]
    # amp2: SpeciesB and SpeciesC (3 SNPs apart) share one cluster and are separated by haplotype_split
    amp2 = summary.loc[(summary.amplicon == "amp2") & summary.consensus.str.contains("_h")]
    assert set(amp2.best_species) == {"SpeciesB", "SpeciesC"} and (amp2.status == "PASS_SPECIES").all()
    assert list(amp2.consensus) == ["Test2_amp2_c1_h1", "Test2_amp2_c1_h2"]
    variants = pd.read_csv(tmp_path / "results/haplotypes/Test2/amp2/variants.tsv", sep="\t")
    assert sorted(variants.position) == [96, 211, 334]
    sizes = pd.read_csv(tmp_path / "results/haplotypes/Test2/amp2/cluster_sizes.tsv", sep="\t").set_index("consensus")
    by_species = dict(zip(amp2.consensus, amp2.best_species))
    reads = {by_species[c]: sizes.loc[c, "reads_mapped"] for c in amp2.consensus}
    assert abs(reads["SpeciesB"] - 50) <= 3 and abs(reads["SpeciesC"] - 30) <= 3
    assert (summary.loc[summary.amplicon == "16S", "best_species"] == "Toyibacter primus").any()
    votes = pd.read_csv(tmp_path / "results/identify/Test2/read_votes.tsv", sep="\t").set_index(["amplicon", "species"])
    assert votes.loc[("amp1", "SpeciesA"), "reads"] == 40 and votes.loc[("amp1", "SpeciesB"), "reads"] == 50
    assert votes.loc[("amp1", "SpeciesD"), "reads"] == 30
    assert votes.loc[("amp2", "SpeciesC"), "reads"] == 30
    novelty = pd.read_csv(tmp_path / "results/identify/Test1/novelty_summary.tsv", sep="\t")
    assert not novelty.is_novel.any()


@pytest.mark.parametrize("header, message", [
    (">my amp|species=SpeciesA", "without spaces"),
    (">amp1 |species=SpeciesA", "without spaces"),
    (">amp1|species=SpeciesA\nACGT\n>amp1|species=SpeciesA", "duplicate sequence ids"),
])
def test_panel_validation(root, tmp_path, header, message):
    panel = tmp_path / "panel.fasta"
    panel.write_text(f"{header}\nACGT\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(f"samples: {root / 'config' / 'samples.tsv'}\npanel:\n  fasta: {panel}\nrrna: {{}}\n")
    r = snakemake(root, tmp_path, "-n", "--cores", "2", "--configfile", str(cfg))
    assert r.returncode != 0 and message in r.stderr + r.stdout
