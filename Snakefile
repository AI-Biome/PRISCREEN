"""
Main Snakefile for ONT Custom Amplicon Analysis Pipeline
"""

import os, re, json, textwrap
import pandas as pd

from snakemake.shell import shell
shell.executable("/bin/bash")

configfile: "config/config.yaml"

# ============================================================================
# Helper Functions
# ============================================================================

def has_flag(f):
    try:
        with open("/proc/cpuinfo") as fh:
            return any(f in line for line in fh)
    except FileNotFoundError:
        return False

def get_amplicons_from_panel(panel_fa):
    amplicons = set()
    with open(panel_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].strip()
                amp = head.split("|")[0]
                amplicons.add(amp)
    return sorted(amplicons)

# ============================================================================
# Configuration & Global Variables
# ============================================================================

# Sample information
SAMPLES_DF = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES = SAMPLES_DF["sample"].tolist()
FASTQ = {r.sample: r.fastq for r in SAMPLES_DF.itertuples()}

# Panel and amplicons
PANEL_FASTA = config["multi_species"]["panel_fasta"]
AMP_LIST = get_amplicons_from_panel(PANEL_FASTA)

# Thread configuration
THREADS_MAP = int(config["threads"].get("map", 8))
THREADS_POLISH = int(config["threads"].get("polish", 8))
THREADS_MMSEQS = int(config["threads"].get("mmseqs", 8))

# Clustering parameters
CLUSTER = bool(config.get("cluster", {}).get("enable", True))
CLUSTER_ID = float(config.get("cluster", {}).get("id_threshold", 0.995))

# Identification parameters
IDENT_MIN_PID = float(config["identify"]["min_pident"])
IDENT_MIN_QCOV = float(config["identify"]["min_qcov"])
IDENT_TOP_DELTA = float(config["identify"]["top_delta"])

# ============================================================================
# Target Rule
# ============================================================================

rule all:
    input:
        expand("results/identify/{sample}/species_summary.tsv", sample=SAMPLES)

# ============================================================================
# Include Rule Modules
# ============================================================================

include: "rules/database.smk"
include: "rules/mapping.smk"
include: "rules/clustering.smk"
include: "rules/consensus.smk"
include: "rules/identification.smk"
