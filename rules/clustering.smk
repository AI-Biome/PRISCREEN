"""
Rules for clustering amplicon reads
"""

rule cluster_vsearch:
    input:
        fq = "results/bin/{sample}/fastq/{amplicon}.fq.gz"
    output:
        cents = "results/cluster/{sample}/{amplicon}/centroids.fasta",
        uc    = "results/cluster/{sample}/{amplicon}/clusters.uc"
    conda:
        "../envs/vsearch.yaml"
    params:
        id = CLUSTER_ID
    resources:
        mem_mb=int(slurm_config['SLURM_ARGS']['mem_of_node']) // int(slurm_config['SLURM_ARGS']['cpus_per_task']),
        runtime=int(slurm_config['SLURM_ARGS']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/cluster/{wildcards.sample}/{wildcards.amplicon}
        zcat {input.fq} \
          | vsearch --cluster_fast - --id {params.id} \
            --centroids {output.cents} --uc {output.uc} --strand both --qmask none --relabel centroid_{wildcards.sample}_{wildcards.amplicon}_
        """
