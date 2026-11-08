#!/usr/bin/env python3
"""Draw the PRISCREEN workflow figure (docs/priscreen_pipeline.png and .svg).

Run from the repository root:  python docs/make_pipeline_figure.py
The figure follows the rules of the workflow; update STEPS when rules change.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

INK, MUTED = "#1a1a1a", "#555555"
COL = {  # (fill, edge) per kind of box
    "in": ("#e9eef7", "#33518a"),
    "prep": ("#eaf1f8", "#2c6aa0"),
    "cons": ("#ebf5ee", "#2e7d4f"),
    "ident": ("#fbf1e4", "#b86e12"),
    "flag": ("#f6e7e7", "#a03a3a"),
    "out": ("#e8f4ec", "#2e7d4f"),
}

# (kind, title, detail, per)   per: what one job covers
STEPS = [
    ("prep", "Map reads to the binning reference", "minimap2 map-ont, 10 secondary hits kept", "per sample"),
    ("prep", "Bin reads by amplicon", "primary alignment decides; samtools region query", "per sample and amplicon"),
    ("cons", "Cluster reads", "vsearch cluster_fast; clusters below min_size dropped", "per sample and amplicon"),
    ("cons", "Polish centroids", "Racon rounds, then Medaka; reads re-assigned and counted", "per sample"),
    ("cons", "Split haplotypes", "reads re-mapped; linked minority SNPs and indels phased; one consensus per haplotype", "per sample"),
    ("ident", "Search consensus sequences", "MMseqs2 vs. panel DB or rRNA DB", "per sample and database"),
    ("ident", "Assign taxa", "identity, coverage and bit-score tie band; lowest common rank", "per amplicon"),
    ("flag", "Flag rare, low-resolution and novel", "read support, assigned rank, identity to closest reference", "per amplicon"),
]
INPUTS = [("ONT reads", "FASTQ per sample"),
          ("PRISTINE panel", "amplicon|species=...  (from PRISTINE)"),
          ("rRNA database", "SILVA / GTDB, lineages (16S, 18S, ...)")]
OUTPUTS = [("read_votes.tsv", "per-read species support"),
           ("species_summary.tsv", "taxon per consensus"),
           ("novelty_summary.tsv", "flags + flagged FASTA")]

W, H = 7.2, 11.9
BX, BW, BH = 0.55, 4.6, 0.72
SIDE_X, SIDE_W = 5.45, 1.6


def box(ax, x, y, w, h, kind, title, detail=None, small=False):
    fc, ec = COL[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03,rounding_size=0.08",
                                facecolor=fc, edgecolor=ec, linewidth=1.4, zorder=2))
    if detail is None:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=8.6 if small else 10.5,
                fontweight="bold", color=INK, zorder=3, wrap=True)
    else:
        ax.text(x + w / 2, y + h * 0.66, title, ha="center", va="center", fontsize=8.6 if small else 10.5,
                fontweight="bold", color=INK, zorder=3)
        ax.text(x + w / 2, y + h * 0.30, detail, ha="center", va="center", fontsize=7.2 if small else 8.2,
                color=MUTED, zorder=3)


def arrow(ax, x0, y0, x1, y1, color="#3a3a3a", style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=12, lw=1.3,
                                 color=color, zorder=4, shrinkA=0, shrinkB=0))


def draw(out_png="docs/priscreen_pipeline.png", out_svg="docs/priscreen_pipeline.svg"):
    fig, ax = plt.subplots(figsize=(W, H), dpi=300)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    ax.text(BX, H - 0.32, "PRISCREEN", fontsize=15, fontweight="bold", color=INK, va="center")
    ax.text(BX + 1.75, H - 0.32, "species identification from ONT amplicon reads", fontsize=9.5, color=MUTED, va="center")

    # inputs row
    y_in = H - 1.35
    iw = (BW + SIDE_W + (SIDE_X - BX - BW)) / 3 - 0.12
    x = BX
    in_centres = []
    for title, detail in INPUTS:
        box(ax, x, y_in, iw, 0.62, "in", title, detail, small=True)
        in_centres.append(x + iw / 2)
        x += iw + 0.18

    # step column
    gap = 0.34
    y = y_in - 0.55 - BH
    centres = []
    for i, (kind, title, detail, per) in enumerate(STEPS, start=1):
        box(ax, BX, y, BW, BH, kind, f"{i}. {title}", detail)
        ax.text(BX + BW + 0.12, y + BH / 2, per, fontsize=7.6, color=MUTED, va="center", style="italic")
        centres.append(y)
        y -= BH + gap
    for i in range(len(centres) - 1):
        arrow(ax, BX + BW / 2, centres[i], BX + BW / 2, centres[i + 1] + BH + 0.01)

    # inputs -> steps
    top = centres[0] + BH
    arrow(ax, in_centres[0], y_in, BX + BW / 2 - 0.6, top + 0.01)                 # reads -> mapping
    arrow(ax, in_centres[1], y_in, BX + BW / 2 + 0.3, top + 0.01, color="#33518a")  # panel -> mapping
    # panel and rRNA DB also feed the search step (dashed, right-hand side)
    y_search = centres[5] + BH * 0.22
    xr = BX + BW + 1.82
    for cx, col in ((in_centres[1], "#33518a"), (in_centres[2], "#33518a")):
        ax.plot([cx, cx], [y_in, y_in - 0.22], color=col, lw=1.0, ls=(0, (3, 2)), zorder=1)
        ax.plot([cx, xr], [y_in - 0.22, y_in - 0.22], color=col, lw=1.0, ls=(0, (3, 2)), zorder=1)
    ax.plot([xr, xr], [y_in - 0.22, y_search], color="#33518a", lw=1.0, ls=(0, (3, 2)), zorder=1)
    arrow(ax, xr, y_search, BX + BW + 0.02, y_search, color="#33518a")
    ax.text(xr + 0.05, (y_in - 0.22 + y_search) / 2 + 0.9, "databases (MMseqs2)", fontsize=7.0, color="#33518a",
            va="center", rotation=90)

    # read votes branch from the mapping step
    y_map = centres[0] + BH / 2
    xv = BX - 0.3
    ax.plot([BX, xv], [y_map, y_map], color="#b86e12", lw=1.0, zorder=1)
    ax.plot([xv, xv], [y_map, 0.55 + 0.31], color="#b86e12", lw=1.0, zorder=1)
    ax.text(xv - 0.05, (y_map + 0.9) / 2, "read votes: primary vs. best other-species alignment",
            fontsize=7.2, color="#b86e12", va="center", ha="right", rotation=90)

    # outputs row
    y_out = 0.55
    x = BX
    out_w = (BW + SIDE_W + (SIDE_X - BX - BW)) / 3 - 0.12
    out_centres = []
    for title, detail in OUTPUTS:
        box(ax, x, y_out, out_w, 0.62, "out", title, detail, small=True)
        out_centres.append(x + out_w / 2)
        x += out_w + 0.18
    bottom = centres[-1]
    arrow(ax, BX + BW / 2 - 0.3, bottom, out_centres[1], y_out + 0.63)
    arrow(ax, BX + BW / 2 + 0.4, bottom, out_centres[2], y_out + 0.63)
    arrow(ax, xv, y_out + 0.31, BX, y_out + 0.31, color="#b86e12")

    fig.savefig(out_png, dpi=300, facecolor="white")
    fig.savefig(out_svg, facecolor="white")
    print(f"wrote {out_png} and {out_svg}")


if __name__ == "__main__":
    draw()
