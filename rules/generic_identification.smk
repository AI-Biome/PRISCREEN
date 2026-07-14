rule generic_mmseqs_search:
    input:
        cons="results/consensus/{sample}/{amplicon}/consensus.fasta",
        db_done="results/mmseqs/genericDB/db.done",
        dbdir=directory("results/mmseqs/genericDB")
    output:
        hits="results/generic_identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    conda:
        "../envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) * THREADS_MMSEQS // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/generic_identify/{wildcards.sample}/{wildcards.amplicon}

        if [ ! -s {input.cons} ]; then
            touch {output.hits}
        else
            mmseqs easy-search \
                {input.cons} \
                {input.dbdir}/genericDB \
                {output.hits} \
                tmp_mmseqs_generic_{wildcards.sample}_{wildcards.amplicon} \
                --search-type 3 \
                --threads {threads} \
                --format-output "query,target,pident,qcov,bits"

            rm -rf tmp_mmseqs_generic_{wildcards.sample}_{wildcards.amplicon}
        fi
        """

rule generic_summarize_hits:
    input:
        hits="results/generic_identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    output:
        tsv="results/generic_identify/{sample}/{amplicon}/summary.tsv"
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

rule generic_sample_summary:
    input:
        summaries=expand("results/generic_identify/{{sample}}/{amplicon}/summary.tsv", amplicon=AMP_LIST)
    output:
        tsv="results/generic_identify/{sample}/generic_summary.tsv"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    script:
        "../scripts/sample_summary.py"
