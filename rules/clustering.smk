"""Clustering of the binned reads of one amplicon into candidate sequence variants."""


rule cluster_reads:
    """vsearch cluster_fast with minimap2-like gap penalties (vsearch's defaults, gap open 20, turn the indels of
    ONT reads into mismatches and push read-to-centroid identities far below the true divergence);
    centroids with fewer than cluster.min_size reads are dropped (not polishable)."""
    input:
        fq="results/bin/{sample}/{amplicon}.fq.gz"
    output:
        cents="results/cluster/{sample}/{amplicon}/centroids.fasta",
        all_cents="results/cluster/{sample}/{amplicon}/centroids_all.fasta",
        uc="results/cluster/{sample}/{amplicon}/clusters.uc"
    log:
        "logs/cluster/{sample}/{amplicon}.log"
    threads: 4
    params:
        id=CLUSTER_ID,
        min_size=CLUSTER_MIN_SIZE
    resources:
        mem_mb=mem_mb(4), runtime=RUNTIME_MIN
    container:
        CONTAINERS["vsearch"]
    shell:
        r"""
        set -euo pipefail
        if [ "$(gzip -dc {input.fq} | head -c 1 | wc -c)" = "0" ]; then
            : > {output.cents}; : > {output.all_cents}; : > {output.uc}
            echo "no reads" > {log}
            exit 0
        fi
        gzip -dc {input.fq} \
          | vsearch --cluster_fast - --id {params.id} --strand both --qmask none --sizeout \
                    --gapopen 4I/4E --gapext 2I/2E \
                    --threads {threads} --relabel {wildcards.sample}_{wildcards.amplicon}_centroid_ \
                    --centroids {output.all_cents} --uc {output.uc} 2> {log}
        if [ -s {output.all_cents} ]; then
            vsearch --sortbysize {output.all_cents} --minsize {params.min_size} --output {output.cents} --quiet 2>> {log}
        else
            : > {output.cents}
        fi
        """
