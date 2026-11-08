"""Optional read filtering, mapping to the reference and binning of reads per amplicon."""


def sample_reads(wc):
    return f"results/filtered/{wc.sample}.fq.gz" if FILTER else FASTQ[wc.sample]


rule filter_reads:
    """Length and mean-quality filter (chopper)."""
    input:
        fq=lambda wc: FASTQ[wc.sample]
    output:
        fq=temp("results/filtered/{sample}.fq.gz")
    log:
        "logs/filter/{sample}.log"
    threads: 4
    params:
        min_len=int(cfg("filter", "min_length", default=200)),
        max_len=lambda wc: f"--maxlength {int(cfg('filter', 'max_length', default=0))}"
                           if int(cfg("filter", "max_length", default=0)) > 0 else "",
        min_q=float(cfg("filter", "min_quality", default=9)),
    resources:
        mem_mb=mem_mb(2), runtime=RUNTIME_MIN
    shell:
        r"""
        set -euo pipefail
        pigz -dcf -p 2 {input.fq} \
          | chopper -q {params.min_q} -l {params.min_len} {params.max_len} --threads {threads} 2> {log} \
          | pigz -p 2 > {output.fq}
        """


rule map_to_reference:
    """Map all reads of a sample to the panel and rRNA bin references; keep secondary hits for read votes."""
    input:
        fq=sample_reads,
        ref="results/panel/mapping_reference.fasta"
    output:
        bam=temp("results/map/{sample}.bam"),
        bai=temp("results/map/{sample}.bam.bai")
    log:
        "logs/map/{sample}.log"
    threads: THREADS_MAP
    params:
        sort_threads=min(4, THREADS_MAP)
    resources:
        mem_mb=mem_mb(THREADS_MAP, 4000), runtime=RUNTIME_MIN
    shell:
        r"""
        set -euo pipefail
        minimap2 -t {threads} -ax map-ont -N 10 --secondary=yes {input.ref} {input.fq} 2> {log} \
          | samtools sort -@ {params.sort_threads} -m 1G -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        """


rule bin_reads:
    """Reads whose primary alignment is to a reference of this amplicon, as FASTQ (one BAM region query, no re-reading)."""
    input:
        bam="results/map/{sample}.bam",
        bai="results/map/{sample}.bam.bai",
        fai="results/panel/mapping_reference.fasta.fai"
    output:
        fq="results/bin/{sample}/{amplicon}.fq.gz",
        bed=temp("results/bin/{sample}/{amplicon}.regions.bed")
    log:
        "logs/bin/{sample}/{amplicon}.log"
    threads: 2
    resources:
        mem_mb=mem_mb(2), runtime=RUNTIME_MIN
    shell:
        r"""
        set -euo pipefail
        awk -F'\t' -v A='{wildcards.amplicon}' 'BEGIN{{OFS="\t"}} index($1, A"|")==1 {{print $1, 0, $2}}' {input.fai} > {output.bed}
        if [ ! -s {output.bed} ]; then
            echo "no reference sequence of amplicon {wildcards.amplicon} in {input.fai}; empty bin" > {log}
            : | pigz > {output.fq}
            exit 0
        fi
        # -F 0x904: no unmapped, secondary or supplementary records; samtools fastq restores read orientation.
        samtools view -u -F 0x904 -M -L {output.bed} {input.bam} 2> {log} \
          | samtools fastq -0 /dev/stdout -s /dev/null - 2>> {log} \
          | pigz -p {threads} > {output.fq}
        """
