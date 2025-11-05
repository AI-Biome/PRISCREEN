"""
Mapping and read binning rules
"""

rule mmseqs_createdb:
    input:
        PANEL_FASTA
    output:
        db = "results/mmseqs/panelDB",
        dbtype = "results/mmseqs/panelDB.dbtype"
    conda:
        "../envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/mmseqs
        mmseqs createdb {input} {output.db}
        """

rule map_to_panel:
    input:
        fq = lambda wc: FASTQ[wc.sample],
        ref = PANEL_FASTA
    output:
        bam = temp("results/map/{sample}.panel.sorted.bam"),
        bai = temp("results/map/{sample}.panel.sorted.bam.bai")
    conda:
        "../envs/minimap2.yaml"
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
        names = "results/bin/{sample}/readlists/{amplicon}.names"
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bin/{wildcards.sample}/readlists

        samtools view {input.bam} \
          | awk -F'\t' '{{print $1"\t"$3}}' \
          | awk -F'\t' 'BEGIN{{OFS="\t"}} $2!="*"{{split($2,a,"|"); print $1,a[1]}}' \
          | awk -v A="{wildcards.amplicon}" '$2==A{{print $1}}' \
          | sort -u > {output.names}
        """

rule amplicon_fastq:
    input:
        fq = lambda wc: FASTQ[wc.sample],
        names = "results/bin/{sample}/readlists/{amplicon}.names"
    output:
        fq = "results/bin/{sample}/fastq/{amplicon}.fq.gz"
    conda:
        "../envs/seqtk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bin/{wildcards.sample}/fastq
        seqtk subseq {input.fq} {input.names} | gzip -c > {output.fq}
        """
