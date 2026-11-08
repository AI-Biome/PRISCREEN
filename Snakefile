"""
PRISCREEN: species identification from ONT amplicon reads.

Reads of every sample are binned by amplicon (PRISTINE custom amplicons and/or
rRNA amplicons such as 16S/18S), clustered, polished into consensus sequences,
split into haplotypes where a cluster mixes close variants, and identified
against the matching database. See README.md.
"""

import os
import sys

import pandas as pd

from snakemake.utils import min_version

min_version("8.0")

__author__ = "Tomas Kudlacek, Katharina J. Hoff"
__version__ = open(os.path.join(workflow.basedir, "VERSION")).read().strip()

configfile: os.path.join(workflow.basedir, "config/config.yaml")

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from pristine_to_panel import amplicon_names, find_amplicon_files  # noqa: E402

include: "rules/common.smk"
include: "rules/panel.smk"
include: "rules/mapping.smk"
include: "rules/clustering.smk"
include: "rules/consensus.smk"
include: "rules/haplotypes.smk"
include: "rules/identification.smk"
include: "rules/novelty.smk"


def final_targets():
    targets = expand(CONS_DIR + "/{sample}/{amplicon}/consensus.fasta", sample=SAMPLES, amplicon=AMP_LIST)
    if IDENTIFY:
        targets += expand("results/identify/{sample}/species_summary.tsv", sample=SAMPLES)
        targets += expand("results/identify/{sample}/novelty_summary.tsv", sample=SAMPLES)
        if PANEL_AMPLICONS:
            targets += expand("results/identify/{sample}/read_votes.tsv", sample=SAMPLES)
    return targets


rule all:
    input:
        final_targets()


onstart:
    print(f"PRISCREEN {__version__}: {len(SAMPLES)} sample(s); panel amplicons: "
          f"{', '.join(PANEL_AMPLICONS) or '-'}; rRNA amplicons: {', '.join(RRNA) or '-'}")
