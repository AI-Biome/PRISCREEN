"""Flags for rare, poorly resolved and potentially novel sequences."""


rule flag_rare_novel:
    input:
        summary="results/identify/{sample}/{amplicon}/summary.tsv",
        cons=CONS_DIR + "/{sample}/{amplicon}/consensus.fasta",
        sizes=CONS_DIR + "/{sample}/{amplicon}/cluster_sizes.tsv"
    output:
        tsv="results/identify/{sample}/{amplicon}/novelty_flags.tsv",
        fasta="results/identify/{sample}/{amplicon}/flagged_sequences.fasta"
    log:
        "logs/flag_rare_novel/{sample}/{amplicon}.log"
    params:
        min_reads=int(cfg("novelty", "min_reads", default=5)),
        min_rel=float(cfg("novelty", "min_rel_abundance", default=0.0)),
        novel_pid=float(cfg("novelty", "novel_pident", default=97.0)),
        novel_qcov=float(cfg("novelty", "novel_qcov", default=0.9)),
        min_resolved_rank=str(cfg("novelty", "min_resolved_rank", default="genus"))
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["python"]
    script:
        "../scripts/flag_rare_novel.py"


rule novelty_summary:
    input:
        expand("results/identify/{{sample}}/{amplicon}/novelty_flags.tsv", amplicon=AMP_LIST)
    output:
        tsv="results/identify/{sample}/novelty_summary.tsv"
    log:
        "logs/novelty_summary/{sample}.log"
    params:
        columns=NOVELTY_COLUMNS
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["python"]
    script:
        "../scripts/sample_summary.py"
