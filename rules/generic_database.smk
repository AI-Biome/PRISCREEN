rule generic_mmseqs_createdb:
    input:
        fasta=config["generic_db"]["fasta"]
    output:
        dbdir=directory("results/mmseqs/genericDB"),
        done="results/mmseqs/genericDB/db.done"
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
        mmseqs createdb {input.fasta} {output.dbdir}/genericDB --dbtype 2
        mmseqs createindex {output.dbdir}/genericDB results/mmseqs/generic_tmp --search-type 3
        touch {output.done}
        """
