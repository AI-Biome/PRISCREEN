"""Haplotype splitting of mixed clusters, one job per sample (all amplicons together)."""


rule haplotype_split:
    """Map the binned reads to the polished consensus sequences, call minority SNPs, phase the reads
    and write one consensus per haplotype (results/haplotypes mirrors results/consensus)."""
    input:
        fqs=expand("results/bin/{{sample}}/{amplicon}.fq.gz", amplicon=AMP_LIST),
        cons=expand("results/consensus/{{sample}}/{amplicon}/consensus.fasta", amplicon=AMP_LIST),
        sizes=expand("results/consensus/{{sample}}/{amplicon}/cluster_sizes.tsv", amplicon=AMP_LIST)
    output:
        cons=expand("results/haplotypes/{{sample}}/{amplicon}/consensus.fasta", amplicon=AMP_LIST),
        sizes=expand("results/haplotypes/{{sample}}/{amplicon}/cluster_sizes.tsv", amplicon=AMP_LIST),
        variants=expand("results/haplotypes/{{sample}}/{amplicon}/variants.tsv", amplicon=AMP_LIST)
    log:
        "logs/haplotypes/{sample}.log"
    threads: THREADS_POLISH
    params:
        amplicons=AMP_LIST,
        workdir="results/haplotypes/{sample}/work",
        phasing={k: cfg("haplotypes", k, default=d) for k, d in HAPLOTYPE_DEFAULTS.items()},
        keep_bam=bool(cfg("haplotypes", "keep_bam", default=False))
    resources:
        mem_mb=mem_mb(THREADS_POLISH, 4000), runtime=RUNTIME_MIN
    container:
        CONTAINERS["medaka"]
    script:
        "../scripts/haplotype_split.py"
