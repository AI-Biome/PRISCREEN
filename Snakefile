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
        mkdir -p results/mmseqs
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
        lists = directory("results/bin/{sample}/readlists")
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
        while read AMP; do
          awk -v A="$AMP" '$2==A{print $1}' {output.lists}/qname_amp.tsv > {output.lists}/"$AMP".names
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
            None
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
    run:
        import os, subprocess, tempfile, gzip
        os.makedirs(f"results/consensus/{wildcards.sample}/{wildcards.amplicon}", exist_ok=True)

        seed_fa = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/seed.fa"
        if params.clustered:

            subprocess.check_call(f"cp {input.cents} {seed_fa}", shell=True)
        else:

            with gzip.open(input.fq, "rt") as fh, open(seed_fa, "w") as out:
                h = fh.readline().strip()
                s = fh.readline().strip()
                fh.readline(); fh.readline()
                if not h.startswith("@"):
                    raise RuntimeError("FASTQ appears empty or malformed for seed.")
                out.write(">seed\n" + s + "\n")


        tmp_sam = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/aln.sam"

        def map_and_sort(template_fa):
            shell(
                "minimap2 -t {threads} -ax map-ont --secondary=no {template_fa} {input.fq} > {tmp_sam}"
            )

        map_and_sort(seed_fa)
        for r in range(params.rounds):
            racon_out = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/racon{r}.fa"

            shell(
                f"{input.racon} -t {threads} {input.fq} {tmp_sam} {seed_fa} > {racon_out}"
            )            

            seed_fa = racon_out
            map_and_sort(seed_fa)

        medaka_dir = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/medaka"
        shell(
            f"rm -rf {medaka_dir}"
        )
        shell(
            f"medaka_consensus -i {input.fq} -d {seed_fa} -o {medaka_dir} -t {threads} -m {params.model}"
        )
        shell(
            f"cp {medaka_dir}/consensus.fasta {output.cons}"
        )

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
        
        mmseqs easy-search {input.cons} {input.db}/panelDB {output.hits} tmp_mmseqs_{wildcards.sample}_{wildcards.amplicon} --threads {threads}
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
    run:
        import pandas as pd
        import os
        cols = ["qseqid","sseqid","pident","qcov","alnlen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bits"]
        df = pd.read_csv(input.hits, sep="\t", header=None, names=cols)
        if df.empty:
            out = pd.DataFrame(columns=["consensus","best_species","pident","qcov","status"])
            out.to_csv(output.tsv, sep="\t", index=False)
            return

        def parse_species(s):
            return s.split("|species=")[-1] if "|species=" in s else "UNK"
        df["species"] = df["sseqid"].map(parse_species)


        df_f = df[(df["pident"] >= params.min_pid*100) & (df["qcov"] >= params.min_qcov*100)].copy()

        out_rows = []
        for q, sub in df_f.groupby("qseqid"):
            sub = sub.sort_values(["pident","qcov","bits"], ascending=[False, False, False]).reset_index(drop=True)
            if len(sub)==0:
                out_rows.append((q, "NO_HIT", None, None, "NO_HIT"))
                continue
            top = sub.iloc[0]
            status = "OK"
            if len(sub) > 1:
                p1 = top["pident"]; p2 = sub.iloc[1]["pident"]
                if p1 > 0 and (p1 - p2)/p1 < params.top_delta:
                    status = "AMBIGUOUS"
            out_rows.append((q, top["species"], float(top["pident"])/100.0, float(top["qcov"])/100.0, status))

        out = pd.DataFrame(out_rows, columns=["consensus","best_species","pident","qcov","status"])
        out.to_csv(output.tsv, sep="\t", index=False)

rule sample_summary:
    input:
        summaries = expand("results/identify/{{sample}}/{amplicon}/summary.tsv", amplicon=AMP_LIST)
    output:
        tsv = "results/identify/{sample}/species_summary.tsv"
    conda:
        "envs/samtools.yaml"
    run:
        import pandas as pd, glob, os
        rows = []
        for p in input.summaries:
            amp = p.split("/")[-2]
            if os.path.exists(p) and os.path.getsize(p)>0:
                df = pd.read_csv(p, sep="\t")
                df.insert(0,"amplicon",amp)
                rows.append(df)
        if rows:
            all_df = pd.concat(rows, ignore_index=True)
        else:
            all_df = pd.DataFrame(columns=["amplicon","consensus","best_species","pident","qcov","status"])
        all_df.to_csv(output.tsv, sep="\t", index=False)
