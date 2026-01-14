rule mmseqs_search:
    input:
        cons="results/consensus/{sample}/{amplicon}/consensus.fasta",
        db_done="results/mmseqs/{amplicon}/panelDB/db.done",
        dbdir=directory("results/mmseqs/{amplicon}/panelDB")
    output:
        hits="results/identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    conda:
        "../envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) * THREADS_MMSEQS // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/identify/{wildcards.sample}/{wildcards.amplicon}
        if [ ! -s {input.cons} ]; then
            touch {output.hits}
        else
            mmseqs easy-search {input.cons} {input.dbdir}/panelDB {output.hits} tmp_mmseqs_{wildcards.sample}_{wildcards.amplicon} --search-type 3 --threads {threads} --format-output "query,target,pident,qcov,bits"
            rm -rf tmp_mmseqs_{wildcards.sample}_{wildcards.amplicon}
        fi
        """

rule summarize_hits:
    input:
        hits="results/identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    output:
        tsv="results/identify/{sample}/{amplicon}/summary.tsv"
    conda:
        "../envs/mmseqs.yaml"
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    params:
        min_pid=IDENT_MIN_PID,
        min_qcov=IDENT_MIN_QCOV,
        top_delta=IDENT_TOP_DELTA
    script:
        "../scripts/summarize_hits.py"

rule sample_summary:
    input:
        summaries=expand("results/identify/{{sample}}/{amplicon}/summary.tsv", amplicon=AMP_LIST)
    output:
        tsv="results/identify/{sample}/species_summary.tsv"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    script:
        "../scripts/sample_summary.py"