"""Consensus polishing of cluster centroids, one job per sample (all amplicons together)."""


rule consensus_polish:
    """Racon rounds then Medaka on the centroids of all amplicons of a sample; one consensus FASTA
    and a cluster_sizes.tsv (cluster size and reads mapped to each consensus) per amplicon."""
    input:
        fqs=expand("results/bin/{{sample}}/{amplicon}.fq.gz", amplicon=AMP_LIST),
        cents=expand("results/cluster/{{sample}}/{amplicon}/centroids.fasta", amplicon=AMP_LIST)
    output:
        cons=expand("results/consensus/{{sample}}/{amplicon}/consensus.fasta", amplicon=AMP_LIST),
        sizes=expand("results/consensus/{{sample}}/{amplicon}/cluster_sizes.tsv", amplicon=AMP_LIST)
    log:
        "logs/consensus/{sample}.log"
    threads: THREADS_POLISH
    params:
        amplicons=AMP_LIST,
        workdir="results/consensus/{sample}/polish",
        rounds=int(cfg("consensus", "racon_rounds", default=2)),
        medaka=bool(cfg("consensus", "medaka", default=True)),
        model=cfg("consensus", "medaka_model") or "",
        bacteria=bool(cfg("consensus", "medaka_bacteria", default=False)),
        keep_intermediate=bool(cfg("consensus", "keep_intermediate", default=False))
    resources:
        mem_mb=mem_mb(THREADS_POLISH, 8000), runtime=RUNTIME_MIN
    script:
        "../scripts/consensus_polish.py"
