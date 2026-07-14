"""
Main Snakefile for ONT Custom Amplicon Analysis Pipeline
"""

import os, re, json, textwrap
import pandas as pd
import configparser

from snakemake.shell import shell
shell.executable("/bin/bash")

# Pipeline version
__author__ = "Tomas Kudlacek, Katharina J. Hoff"
__version__ = "0.1.0-beta.1"

configfile: "config/config.yaml"

# Read SLURM configuration from config.ini
slurm_config = configparser.ConfigParser()
slurm_config.read('config.ini')

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

# Taxonomy mode
TAXONOMY_MODE = config.get("taxonomy", {}).get("mode", "panel")

if TAXONOMY_MODE not in {"panel", "hybrid", "none"}:
    raise ValueError("taxonomy.mode must be one of: panel, hybrid, none")

# Panel and amplicons
AMPLICONS = config["amplicons"]
AMP_LIST = sorted(AMPLICONS.keys())
PANEL_FASTA = {a: AMPLICONS[a]["panel_fasta"] for a in AMP_LIST}

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

IDENTIFY = bool(config.get("identify", {}).get("enable", True))
IDENTIFY_DIR = "results/identify" if TAXONOMY_MODE == "panel" else "results/generic_identify"

# ============================================================================
# Target Rule
# ============================================================================

if TAXONOMY_MODE == "panel":
    rule all:
        input:
            expand("results/identify/{sample}/species_summary.tsv", sample=SAMPLES),
            expand("results/identify/{sample}/novelty_summary.tsv", sample=SAMPLES)

elif TAXONOMY_MODE == "hybrid":
    rule all:
        input:
            expand("results/generic_identify/{sample}/generic_summary.tsv", sample=SAMPLES),
            expand("results/generic_identify/{sample}/novelty_summary.tsv", sample=SAMPLES)

else:
    rule all:
        input:
            expand(
                "results/consensus/{sample}/{amplicon}/consensus.fasta",
                sample=SAMPLES,
                amplicon=AMP_LIST
            )

# ============================================================================
# Include Rule Modules
# ============================================================================

include: "rules/mapping.smk"
include: "rules/clustering.smk"
include: "rules/consensus.smk"

if TAXONOMY_MODE == "panel":
    include: "rules/database.smk"
    include: "rules/identification.smk"
    include: "rules/novelty.smk"

elif TAXONOMY_MODE == "hybrid":
    include: "rules/generic_database.smk"
    include: "rules/generic_identification.smk"
    include: "rules/novelty.smk"

# ============================================================================
# Cleanup Handlers
# ============================================================================

onsuccess:
    shell("rm -rf tmp_mmseqs_* 2>/dev/null || true")

onerror:
    shell("rm -rf tmp_mmseqs_* 2>/dev/null || true")
