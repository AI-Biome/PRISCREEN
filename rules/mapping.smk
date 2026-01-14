rule map_to_panel:
    input:
        fq=lambda wc: FASTQ[wc.sample],
        ref=lambda wc: PANEL_FASTA[wc.amplicon]
    output:
        bam=temp("results/map/{sample}/{amplicon}.panel.sorted.bam"),
        bai=temp("results/map/{sample}/{amplicon}.panel.sorted.bam.bai")
    conda:
        "../envs/minimap2.yaml"
    threads: THREADS_MAP
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) * THREADS_MAP // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/map/{wildcards.sample}
        minimap2 -t {threads} -ax map-ont -N 10 --secondary=yes {input.ref} {input.fq} \
          | samtools view -b - \
          | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule amplicon_readnames:
    input:
        bam="results/map/{sample}/{amplicon}.panel.sorted.bam"
    output:
        names="results/bin/{sample}/readlists/{amplicon}.names"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bin/{wildcards.sample}/readlists
        samtools view {input.bam} \
          | awk -F'\t' '$3!="*"{{print $1}}' \
          | sort -u > {output.names}
        """

rule amplicon_fastq:
    input:
        fq=lambda wc: FASTQ[wc.sample],
        names="results/bin/{sample}/readlists/{amplicon}.names"
    output:
        fq="results/bin/{sample}/fastq/{amplicon}.fq.gz"
    conda:
        "../envs/seqtk.yaml"
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
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
