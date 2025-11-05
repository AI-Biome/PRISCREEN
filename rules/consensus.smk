"""
Consensus generation and polishing rules
"""

rule consensus_polish:
    input:
        fq = "results/bin/{sample}/fastq/{amplicon}.fq.gz",
        cents = lambda wc: (
            f"results/cluster/{wc.sample}/{wc.amplicon}/centroids.fasta"
            if CLUSTER else
            None
        )
    output:
        cons = "results/consensus/{sample}/{amplicon}/consensus.fasta"
    conda:
        "../envs/consensus.yaml"
    threads: THREADS_POLISH
    params:
        rounds = int(config["consensus"]["racon_rounds"]),
        model  = config["consensus"]["medaka_model"],
        clustered = CLUSTER
    script:
        "../scripts/consensus_polish.py"
