"""Similarity search of consensus sequences and taxonomic assignment."""

SUMMARY_COLUMNS = ["consensus", "best_species", "assigned_rank", "pident", "qcov", "bits",
                   "n_hits_considered", "n_distinct_species", "closest_target", "status"]
NOVELTY_COLUMNS = SUMMARY_COLUMNS + ["read_count", "cluster_size", "separated_by", "rel_abundance", "is_rare",
                                     "is_low_resolution", "is_novel", "flag_reason", "blast_url"]


def amplicons_of_db(db):
    return PANEL_AMPLICONS if db == "panel" else [db]


rule mmseqs_search:
    """Search the consensus sequences of all amplicons that use one database (panel or an rRNA database)."""
    input:
        cons=lambda wc: expand(CONS_DIR + "/{sample}/{amplicon}/consensus.fasta",
                               sample=wc.sample, amplicon=amplicons_of_db(wc.db)),
        db_done="results/db/{db}/db.done"
    output:
        query="results/search/{sample}/{db}/query.fasta",
        hits="results/search/{sample}/{db}/mmseqs_hits.tsv",
        tmp=temp(directory("results/search/{sample}/{db}/tmp_mmseqs"))
    log:
        "logs/search/{sample}/{db}.log"
    threads: THREADS_MMSEQS
    params:
        db="results/db/{db}/db",
        max_seqs=int(cfg("identify", "max_seqs", default=300))
    resources:
        mem_mb=mem_mb(THREADS_MMSEQS, 4000), runtime=RUNTIME_MIN
    container:
        CONTAINERS["mmseqs2"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.tmp}
        cat {input.cons} > {output.query}
        if [ ! -s {output.query} ]; then
            : > {output.hits}
            echo "no consensus sequences" > {log}
            exit 0
        fi
        mmseqs easy-search {output.query} {params.db} {output.hits} {output.tmp} \
            --search-type 3 --threads {threads} --max-seqs {params.max_seqs} -v 1 \
            --format-output "query,target,theader,pident,qcov,tcov,alnlen,bits,evalue" > {log} 2>&1
        """


rule summarize_hits:
    """Taxon per consensus of one amplicon from the hits of its database (consensus without hits: NO_HIT)."""
    input:
        hits=lambda wc: f"results/search/{wc.sample}/{database_of(wc.amplicon)}/mmseqs_hits.tsv",
        cons=CONS_DIR + "/{sample}/{amplicon}/consensus.fasta"
    output:
        tsv="results/identify/{sample}/{amplicon}/summary.tsv"
    log:
        "logs/summarize_hits/{sample}/{amplicon}.log"
    params:
        min_pident=float(cfg("identify", "min_pident", default=98.5)),
        min_qcov=float(cfg("identify", "min_qcov", default=0.9)),
        top_delta=float(cfg("identify", "top_delta", default=0.01))
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["python"]
    script:
        "../scripts/summarize_hits.py"


rule species_summary:
    input:
        expand("results/identify/{{sample}}/{amplicon}/summary.tsv", amplicon=AMP_LIST)
    output:
        tsv="results/identify/{sample}/species_summary.tsv"
    log:
        "logs/species_summary/{sample}.log"
    params:
        columns=SUMMARY_COLUMNS
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["python"]
    script:
        "../scripts/sample_summary.py"


rule read_votes:
    """Per-read species support for panel amplicons from the binning alignment (primary vs. best other-species hit)."""
    input:
        bam="results/map/{sample}.bam",
        bai="results/map/{sample}.bam.bai"
    output:
        tsv="results/identify/{sample}/read_votes.tsv"
    log:
        "logs/read_votes/{sample}.log"
    params:
        min_margin=int(cfg("identify", "vote_min_margin", default=4)),
        panel_amplicons=PANEL_AMPLICONS
    resources:
        mem_mb=mem_mb(1, 2000), runtime=RUNTIME_MIN
    container:
        CONTAINERS["medaka"]
    script:
        "../scripts/read_votes.py"
