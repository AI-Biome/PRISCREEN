"""
Rules for database preparation
"""

rule mmseqs_createdb:
    input:
        PANEL_FASTA
    output:
        db = directory("results/mmseqs/panelDB")
    conda:
        "../envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.db}
        mmseqs createdb {input} {output.db}/panelDB
        """
