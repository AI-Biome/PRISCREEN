rule flag_rare_novel:
    input:
        summary="results/identify/{sample}/{amplicon}/summary.tsv",
        cons="results/consensus/{sample}/{amplicon}/consensus.fasta",
        uc=lambda wc: f"results/cluster/{wc.sample}/{wc.amplicon}/clusters.uc" if CLUSTER else ""
    output:
        tsv="results/identify/{sample}/{amplicon}/novelty_flags.tsv",
        fasta="results/identify/{sample}/{amplicon}/flagged_sequences.fasta"
    conda:
        "../envs/mmseqs.yaml"
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    params:
        min_reads=int(config.get("novelty", {}).get("min_reads", 5)),
        min_rel=float(config.get("novelty", {}).get("min_rel_abundance", 0.0)),
        novel_pid=float(config.get("novelty", {}).get("novel_pident", 97.0)),
        novel_qcov=float(config.get("novelty", {}).get("novel_qcov", 0.90)),
        treat_unclassified_as_novel=bool(config.get("novelty", {}).get("treat_unclassified_as_novel", True)),
        min_resolved_rank=str(config.get("novelty", {}).get("min_resolved_rank", "genus"))
    script:
        "../scripts/flag_rare_novel.py"


rule sample_novelty_summary:
    input:
        flags=expand("results/identify/{{sample}}/{amplicon}/novelty_flags.tsv", amplicon=AMP_LIST)
    output:
        tsv="results/identify/{sample}/novelty_summary.tsv"
    conda:
        "../envs/mmseqs.yaml"
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    script:
        "../scripts/sample_novelty_summary.py"