"""Unit tests of the Snakemake-run scripts (pure functions, no tools needed)."""
import pandas as pd

import summarize_hits
from flag_rare_novel import flag
from read_votes import parse_sam as _parse_sam


def parse_sam(lines, min_margin, panel_amplicons=None):
    """Vote rows as a DataFrame indexed by (amplicon, species), for compact assertions."""
    rows = _parse_sam(lines, min_margin, panel_amplicons)
    cols = ["amplicon", "species", "reads", "reads_unique", f"reads_margin_ge_{min_margin}", "fraction_of_amplicon_reads"]
    return pd.DataFrame(rows, columns=cols)
from taxonomy import lowest_common_rank, parse_header_taxonomy, rank_index


def hits_df(rows):
    df = pd.DataFrame(rows, columns=summarize_hits.COLUMNS)
    for c in ["pident", "qcov", "tcov", "alnlen", "bits", "evalue"]:
        df[c] = pd.to_numeric(df[c])
    return df


def summarize(rows, queries, min_pident=98.5, min_qcov=0.9, top_delta=0.01):
    return summarize_hits.summarize(hits_df(rows), queries, min_pident, min_qcov, top_delta,
                                    parse_header_taxonomy, lowest_common_rank)


def test_species_call_and_tie_band():
    rows = [["c1", "amp1|species=SpA", "amp1|species=SpA", 100.0, 1.0, 1.0, 400, 718, 1e-100],
            ["c1", "amp1|species=SpB", "amp1|species=SpB", 99.0, 1.0, 1.0, 400, 700, 1e-90]]
    out = summarize(rows, ["c1"])
    r = out.iloc[0]
    assert r.status == "PASS_SPECIES" and r.best_species == "SpA" and r.n_hits_considered == 1
    # with a wide tie band both species are tied: no common rank for single-token names
    r = summarize(rows, ["c1"], top_delta=0.1).iloc[0]
    assert r.status == "AMBIGUOUS" and r.n_distinct_species == 2
    # a hit below min_pident never enters the tie band
    rows[1][3] = 98.0
    assert summarize(rows, ["c1"], top_delta=0.1).iloc[0].status == "PASS_SPECIES"


def test_genus_level_when_species_tie():
    rows = [["c1", "l|species=Yersinia_pestis", "l|species=Yersinia_pestis", 99.9, 1.0, 1.0, 1000, 1800, 0],
            ["c1", "l|species=Yersinia_pseudotuberculosis", "l|species=Yersinia_pseudotuberculosis", 99.9, 1.0, 1.0, 1000, 1800, 0]]
    r = summarize(rows, ["c1"]).iloc[0]
    assert r.status == "PASS_GENUS" and r.best_species == "Yersinia" and r.assigned_rank == "genus"


def test_below_threshold_and_no_hit_rows():
    rows = [["c1", "t", "t Bacteria;P;C;O;F;G;G sp", 96.0, 1.0, 1.0, 400, 600, 0]]
    out = summarize(rows, ["c1", "c2"])
    assert list(out.status) == ["BELOW_THRESHOLD", "NO_HIT"]
    assert out.iloc[0].closest_target == "t" and out.iloc[0].pident == 96.0
    assert out.iloc[1].pident == ""


def test_qcov_threshold():
    rows = [["c1", "t", "t|species=X", 100.0, 0.5, 1.0, 200, 300, 0]]
    assert summarize(rows, ["c1"]).iloc[0].status == "BELOW_THRESHOLD"


def test_flags():
    df = pd.DataFrame([
        dict(consensus="c1", best_species="SpA", assigned_rank="species", pident=100.0, qcov=1.0, bits=700, status="PASS_SPECIES"),
        dict(consensus="c2", best_species="Fam", assigned_rank="family", pident=99.0, qcov=1.0, bits=600, status="PASS_FAMILY"),
        dict(consensus="c3", best_species="", assigned_rank="", pident=95.0, qcov=1.0, bits=500, status="BELOW_THRESHOLD"),
        dict(consensus="c4", best_species="", assigned_rank="", pident="", qcov="", bits="", status="NO_HIT"),
    ])
    sizes = {"c1": (40, 50, ""), "c2": (10, 12, "indel"), "c3": (3, 2, "snp"), "c4": (3, 0, "")}
    out = flag(df, sizes, min_reads=5, min_rel=0.0, novel_pid=97.0, novel_qcov=0.9,
               min_resolved_rank="genus", rank_index=rank_index)
    out = out.set_index("consensus")
    assert out.loc["c1", "read_count"] == 50 and out.loc["c1", "cluster_size"] == 40
    assert not out.loc["c1", ["is_rare", "is_low_resolution", "is_novel"]].any()
    assert out.loc["c2", "is_low_resolution"] and out.loc["c2", "flag_reason"] == "low_resolution<genus;indel_only_split"
    assert not out.loc["c2", "is_novel"] and out.loc["c2", "separated_by"] == "indel"
    assert out.loc["c3", "is_rare"] and out.loc["c3", "is_novel"]
    assert "below_threshold" in out.loc["c3", "flag_reason"] and "low_pid<97.0" in out.loc["c3", "flag_reason"]
    assert out.loc["c4", "is_novel"] and out.loc["c4", "is_rare"] and out.loc["c4", "flag_reason"] == "rare;no_hit"


def sam_line(q, flag, ref, score):
    return "\t".join([q, str(flag), ref, "1", "60", "400M", "*", "0", "0", "*", "*", f"AS:i:{score}"])


def test_read_votes():
    lines = [
        sam_line("r1", 0, "amp1|species=SpA", 700), sam_line("r1", 256, "amp1|species=SpB", 690),
        sam_line("r2", 16, "amp1|species=SpA", 700), sam_line("r2", 272, "amp1|species=SpB", 698),  # margin 2
        sam_line("r3", 0, "amp1|species=SpB", 650), sam_line("r3", 256, "amp1|species=SpA", 650),   # tie
        sam_line("r4", 0, "amp2|species=SpA", 600), sam_line("r4", 256, "amp1|species=SpA", 300),   # other amplicon
        sam_line("r5", 4, "*", 0),
        sam_line("r6", 0, "16S|TOY0001", 900),
        sam_line("r7", 2048, "amp1|species=SpA", 100),
    ]
    df = parse_sam(lines, min_margin=4).set_index(["amplicon", "species"])
    assert df.loc[("amp1", "SpA"), "reads"] == 2
    assert df.loc[("amp1", "SpA"), "reads_unique"] == 2
    assert df.loc[("amp1", "SpA"), "reads_margin_ge_4"] == 1
    assert df.loc[("amp1", "SpB"), "reads"] == 1 and df.loc[("amp1", "SpB"), "reads_unique"] == 0
    assert df.loc[("amp2", "SpA"), "reads_unique"] == 1
    assert ("16S", None) not in df.index and len(df) == 3


def test_read_votes_order_independent_and_panel_only():
    """Coordinate-sorted BAMs list a secondary record before the primary of the same read."""
    lines = [
        sam_line("r1", 256, "amp1|species=SpA|genome=1", 700),   # secondary of the same species, seen first
        sam_line("r1", 0, "amp1|species=SpA|genome=2", 700),
        sam_line("r2", 272, "amp1|species=SpB", 698),            # other species, seen first: margin 2
        sam_line("r2", 16, "amp1|species=SpA|genome=2", 700),
        sam_line("r3", 0, "16S|species=Foo", 500),               # rRNA bin reference with a species= field
    ]
    df = parse_sam(lines, min_margin=4, panel_amplicons={"amp1"}).set_index(["amplicon", "species"])
    assert df.loc[("amp1", "SpA"), "reads"] == 2
    assert df.loc[("amp1", "SpA"), "reads_unique"] == 2 and df.loc[("amp1", "SpA"), "reads_margin_ge_4"] == 1
    assert ("16S", "Foo") not in df.index
    assert ("16S", "Foo") in parse_sam(lines, min_margin=4).set_index(["amplicon", "species"]).index
