"""
Rules for species identification
"""

rule mmseqs_search:
    input:
        cons = "results/consensus/{sample}/{amplicon}/consensus.fasta",
        db   = rules.mmseqs_createdb.output.db
    output:
        hits = "results/identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    conda:
        "../envs/mmseqs.yaml"
    threads: THREADS_MMSEQS
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/identify/{wildcards.sample}/{wildcards.amplicon}

        # Check if consensus file is empty, if so create empty output
        if [ ! -s {input.cons} ]; then
            touch {output.hits}
        else
            mmseqs easy-search {input.cons} {input.db}/panelDB {output.hits} tmp_mmseqs_{wildcards.sample}_{wildcards.amplicon} --search-type 3 --threads {threads} --format-output "query,target,pident,qcov,bits"
        fi
        """

rule summarize_hits:
    input:
        hits = "results/identify/{sample}/{amplicon}/mmseqs_hits.tsv"
    output:
        tsv = "results/identify/{sample}/{amplicon}/summary.tsv"
    conda:
        "../envs/mmseqs.yaml"
    params:
        min_pid = IDENT_MIN_PID,
        min_qcov = IDENT_MIN_QCOV,
        top_delta = IDENT_TOP_DELTA
    script:
        "../scripts/summarize_hits.py"

rule sample_summary:
    input:
        summaries = expand("results/identify/{{sample}}/{amplicon}/summary.tsv", amplicon=AMP_LIST)
    output:
        tsv = "results/identify/{sample}/species_summary.tsv"
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/sample_summary.py"
