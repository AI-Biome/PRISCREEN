"""
Rules for database preparation
"""

rule mmseqs_createdb:
    input:
        panel=lambda wc: PANEL_FASTA[wc.amplicon]
    output:
        dbdir=directory("results/mmseqs/{amplicon}/panelDB"),
        done="results/mmseqs/{amplicon}/panelDB/db.done"
    conda:
        "../envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) * THREADS_MMSEQS // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.dbdir}
        mmseqs createdb {input.panel} {output.dbdir}/panelDB
        mmseqs createindex {output.dbdir}/panelDB results/mmseqs/tmp --search-type 3
        touch {output.done}
        """
