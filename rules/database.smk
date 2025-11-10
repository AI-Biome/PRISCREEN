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
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.db}
        mmseqs createdb {input} {output.db}/panelDB
        """
