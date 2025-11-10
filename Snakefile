import os, re, json, textwrap
import pandas as pd

from snakemake.shell import shell
shell.executable("/bin/bash")

configfile: "config/config.yaml"

def has_flag(f):
    try:
        with open("/proc/cpuinfo") as fh:
            return any(f in line for line in fh)
    except FileNotFoundError:
        return False

# if has_flag("avx2"):
#     RACON_BIN = "software/racon/bin/racon_avx2"
# else:
#     RACON_BIN = "software/racon/bin/racon_generic"

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

rule mmseqs_createdb:
    input:
        PANEL_FASTA
    output:
        db = directory("results/mmseqs/panelDB")
    conda:
        "envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.db}
        mmseqs createdb {input} {output.db}/panelDB
        """

rule map_to_panel:
    input:
        fq = lambda wc: FASTQ[wc.sample],
        ref = PANEL_FASTA
    output:
        bam = temp("results/map/{sample}.panel.sorted.bam"),
        bai = temp("results/map/{sample}.panel.sorted.bam.bai")
    conda:
        "envs/minimap2.yaml"
    threads: THREADS_MAP
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/map
        minimap2 -t {threads} -ax map-ont -N 10 --secondary=yes {input.ref} {input.fq} \
          | samtools view -b - \
          | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule amplicon_readnames:
    input:
        bam = "results/map/{sample}.panel.sorted.bam"
    output:
        lists = directory("results/bin/{sample}/readlists"),
        names = expand("results/bin/{{sample}}/readlists/{amp}.names", amp=AMP_LIST)
    params:
        amps = " ".join(AMP_LIST)
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.lists}

        samtools view {input.bam} \
          | awk -F'\t' '{{print $1"\t"$3}}' \
          | awk -F'\t' 'BEGIN{{OFS="\t"}} $2!="*"{{split($2,a,"|"); print $1,a[1]}}' \
          | sort -u > {output.lists}/qname_amp.tsv

        cut -f2 {output.lists}/qname_amp.tsv | sort -u > {output.lists}/amps.txt

        # Create empty .names files for all expected amplicons
        for AMP in {params.amps}; do
          touch {output.lists}/"$AMP".names
        done

        # Fill in the .names files for amplicons that have reads
        while read AMP; do
          awk -v A="$AMP" '$2==A{{print $1}}' {output.lists}/qname_amp.tsv > {output.lists}/"$AMP".names
        done < {output.lists}/amps.txt
        """

rule amplicon_fastq:
    input:
        fq = lambda wc: FASTQ[wc.sample],
        names = "results/bin/{sample}/readlists/{amplicon}.names"
    output:
        fq = "results/bin/{sample}/fastq/{amplicon}.fq.gz"
    conda:
        "envs/seqtk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bin/{wildcards.sample}/fastq
        if [ -s {input.names} ]; then
          seqtk subseq {input.fq} {input.names} | gzip -c > {output.fq}
        else
          printf "" | gzip -c > {output.fq}
        fi
        """

rule cluster_vsearch:
    input:
        fq = "results/bin/{sample}/fastq/{amplicon}.fq.gz"
    output:
        cents = "results/cluster/{sample}/{amplicon}/centroids.fasta",
        uc    = "results/cluster/{sample}/{amplicon}/clusters.uc"
    conda:
        "envs/vsearch.yaml"
    params:
        id = CLUSTER_ID
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/cluster/{wildcards.sample}/{wildcards.amplicon}
        zcat {input.fq} \
          | vsearch --cluster_fast - --id {params.id} \
            --centroids {output.cents} --uc {output.uc} --strand both --qmask none --relabel centroid_{wildcards.sample}_{wildcards.amplicon}_
        """

rule build_racon:
    output:
        generic = "software/racon_build/generic/bin/racon",
        avx2    = "software/racon_build/avx2/bin/racon"
    conda:
        "envs/racon-build.yaml"
    threads: 4
    shell:
        r"""
        set -euo pipefail
        rm -rf software/racon_src software/racon_build
        git clone --recursive https://github.com/lbcb-sci/racon.git software/racon_src

        # Generic (portable)
        mkdir -p software/racon_build/generic
        cd software/racon_build/generic
        cmake -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_CXX_FLAGS="-O3 -march=x86-64 -mtune=generic" \
              -DCMAKE_C_FLAGS="-O3 -march=x86-64 -mtune=generic" \
              ../../racon_src
        make -j {threads}
        install -D bin/racon ../../../racon/bin/racon_generic

        # AVX2 (fast)
        cd ..
        mkdir -p avx2
        cd avx2
        cmake -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_CXX_FLAGS="-O3 -mavx2 -mtune=native" \
              -DCMAKE_C_FLAGS="-O3 -mavx2 -mtune=native" \
              ../../racon_src
        make -j {threads}
        install -D bin/racon ../../../racon/bin/racon_avx2
        """

def select_racon_output():
    outs = rules.build_racon.output
    if has_flag("avx2"):
        return outs.avx2

    return outs.generic

rule consensus_polish:
    input:
        racon = select_racon_output(),
        fq = "results/bin/{sample}/fastq/{amplicon}.fq.gz",
        cents = lambda wc: (
            f"results/cluster/{wc.sample}/{wc.amplicon}/centroids.fasta"
            if CLUSTER else
            []
        )
    output:
        cons = "results/consensus/{sample}/{amplicon}/consensus.fasta"
    conda:
        "envs/consensus.yaml"
    threads: THREADS_POLISH
    params:
        rounds = int(config["consensus"]["racon_rounds"]),
        model  = config["consensus"]["medaka_model"],
        clustered = CLUSTER
    script:
        "scripts/consensus_polish.py"

rule mmseqs_search:
    input:
        cons = "results/consensus/{sample}/{amplicon}/consensus.fasta",
        db   = rules.mmseqs_createdb.output.db
    output:
        hits = "results/identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    conda:
        "envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/identify/{wildcards.sample}/{wildcards.amplicon}

        # Check if consensus file is empty, if so create empty output
        if [ ! -s {input.cons} ]; then
            touch {output.hits}
        else
            mmseqs easy-search {input.cons} {input.db}/panelDB {output.hits} tmp_mmseqs_{wildcards.sample}_{wildcards.amplicon} --search-type 3 --threads {threads}
        fi
        """

rule summarize_hits:
    input:
        hits = "results/identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    output:
        tsv = "results/identify/{sample}/{amplicon}/summary.tsv"
    conda:
        "envs/mmseqs.yaml"
    params:
        min_pid = IDENT_MIN_PID,
        min_qcov = IDENT_MIN_QCOV,
        top_delta = IDENT_TOP_DELTA
    script:
        "scripts/summarize_hits.py"

rule sample_summary:
    input:
        summaries = expand("results/identify/{{sample}}/{amplicon}/summary.tsv", amplicon=AMP_LIST)
    output:
        tsv = "results/identify/{sample}/species_summary.tsv"
    conda:
        "envs/samtools.yaml"
    script:
        "scripts/sample_summary.py"
