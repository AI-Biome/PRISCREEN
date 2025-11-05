import os, re, json, textwrap
import pandas as pd

configfile: "config/config.yaml"

SAMPLES_DF = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES = SAMPLES_DF["sample"].tolist()
FASTQ = {r.sample: r.fastq for r in SAMPLES_DF.itertuples()}

PANEL_FASTA = config["multi_species"]["panel_fasta"]

THREADS_MAP = int(config["threads"].get("map", 8))
THREADS_POLISH = int(config["threads"].get("polish", 8))
THREADS_MMSEQS = int(config["threads"].get("mmseqs", 8))

def get_amplicons_from_panel(panel_fa):
    amplicons = set()
    with open(panel_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].strip()
                amp = head.split("|")[0]
                amplicons.add(amp)
    return sorted(amplicons)

AMP_LIST = get_amplicons_from_panel(PANEL_FASTA)

CLUSTER = bool(config.get("cluster", {}).get("enable", True))
CLUSTER_ID = float(config.get("cluster", {}).get("id_threshold", 0.995))

IDENT_MIN_PID = float(config["identify"]["min_pident"])
IDENT_MIN_QCOV = float(config["identify"]["min_qcov"])
IDENT_TOP_DELTA = float(config["identify"]["top_delta"])

rule all:
    input:
        expand("results/identify/{sample}/species_summary.tsv", sample=SAMPLES)

# Include modular rule files
include: "rules/mapping.smk"
include: "rules/clustering.smk"
include: "rules/consensus.smk"
include: "rules/identification.smk"
