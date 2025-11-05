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
    run:
        import os, subprocess, tempfile, gzip
        os.makedirs(f"results/consensus/{wildcards.sample}/{wildcards.amplicon}", exist_ok=True)

        seed_fa = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/seed.fa"
        if params.clustered:

            subprocess.check_call(f"cp {input.cents} {seed_fa}", shell=True)
        else:

            with gzip.open(input.fq, "rt") as fh, open(seed_fa, "w") as out:
                h = fh.readline().strip()
                s = fh.readline().strip()
                fh.readline(); fh.readline()
                if not h.startswith("@"):
                    raise RuntimeError("FASTQ appears empty or malformed for seed.")
                out.write(">seed\n" + s + "\n")


        tmp_bam = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/aln.bam"

        def map_and_sort(template_fa):
            cmd = f"minimap2 -t {threads} -ax map-ont {template_fa} {input.fq} | samtools view -b - | samtools sort -o {tmp_bam}"
            subprocess.check_call(cmd, shell=True)
            subprocess.check_call(f"samtools index {tmp_bam}", shell=True)

        map_and_sort(seed_fa)
        for r in range(params.rounds):
            racon_out = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/racon{r}.fa"

            cmd = f"racon -t {threads} <(zcat {input.fq}) {tmp_bam} {seed_fa} > {racon_out}"
            subprocess.check_call(["bash","-lc",cmd])
            seed_fa = racon_out
            map_and_sort(seed_fa)


        medaka_dir = f"results/consensus/{wildcards.sample}/{wildcards.amplicon}/medaka"
        subprocess.check_call(f"rm -rf {medaka_dir}", shell=True)
        cmd = f"medaka_consensus -i {input.fq} -d {seed_fa} -o {medaka_dir} -t {threads} -m {params.model}"
        subprocess.check_call(cmd, shell=True)

        subprocess.check_call(f"cp {medaka_dir}/consensus.fasta {output.cons}", shell=True)
