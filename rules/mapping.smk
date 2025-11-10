"""
Rules for mapping reads to panel and binning by amplicon
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
        lists = directory("results/bin/{sample}/readlists"),
        names = expand("results/bin/{{sample}}/readlists/{amp}.names", amp=AMP_LIST)
    params:
        amps = " ".join(AMP_LIST)
    conda:
        "../envs/samtools.yaml"
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
        "../envs/seqtk.yaml"
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
