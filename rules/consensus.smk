"""Consensus polishing of cluster centroids, per sample: Racon rounds, then Medaka.

Every tool runs in its own public container, so the step is a chain of small rules per
sample (the seeds of all amplicons are polished together, Medaka runs once per sample):
polish_prepare -> racon_map / racon (x racon_rounds) -> medaka (or map_to_draft) ->
count_mapped -> consensus_collect.
"""

RACON_ROUNDS = int(cfg("consensus", "racon_rounds", default=2))
USE_MEDAKA = bool(cfg("consensus", "medaka", default=True))
POLISH = "results/consensus/{sample}/polish"


wildcard_constraints:
    round="|".join(str(r) for r in range(1, RACON_ROUNDS + 1)) or "none",


rule polish_prepare:
    """Seeds of all amplicons of a sample (centroids with >= min_size reads, renamed <sample>_<amplicon>_c<i>)
    and one FASTQ with the reads of all bins."""
    input:
        fqs=expand("results/bin/{{sample}}/{amplicon}.fq.gz", amplicon=AMP_LIST),
        cents=expand("results/cluster/{{sample}}/{amplicon}/centroids.fasta", amplicon=AMP_LIST)
    output:
        seeds=POLISH + "/seeds.fasta",
        table=POLISH + "/seeds.tsv",
        reads=temp(POLISH + "/reads.fq.gz")
    log:
        "logs/consensus/{sample}_prepare.log"
    params:
        amplicons=AMP_LIST
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["python"]
    script:
        "../scripts/prepare_polish.py"


def racon_draft(wc):
    r = int(wc.round)
    return POLISH.format(sample=wc.sample) + ("/seeds.fasta" if r == 1 else f"/racon{r - 1}.fasta")


rule racon_map:
    """Reads mapped to the current draft (primary alignments only, so each read polishes one seed)."""
    input:
        draft=racon_draft,
        reads=POLISH + "/reads.fq.gz"
    output:
        sam=temp(POLISH + "/racon{round}.sam")
    log:
        "logs/consensus/{sample}_racon{round}_map.log"
    threads: THREADS_POLISH
    resources:
        mem_mb=mem_mb(THREADS_POLISH, 4000), runtime=RUNTIME_MIN
    container:
        CONTAINERS["minimap2"]
    shell:
        r"""
        set -euo pipefail
        if [ -s {input.draft} ] && [ -s {input.reads} ]; then
            minimap2 -t {threads} -ax map-ont --secondary=no {input.draft} {input.reads} > {output.sam} 2> {log}
        else
            : > {output.sam}; echo "no seeds or no reads" > {log}
        fi
        """


rule racon:
    input:
        draft=racon_draft,
        sam=POLISH + "/racon{round}.sam",
        reads=POLISH + "/reads.fq.gz"
    output:
        fa=POLISH + "/racon{round}.fasta"
    log:
        "logs/consensus/{sample}_racon{round}.log"
    threads: THREADS_POLISH
    resources:
        mem_mb=mem_mb(THREADS_POLISH, 4000), runtime=RUNTIME_MIN
    container:
        CONTAINERS["racon"]
    shell:
        r"""
        set -euo pipefail
        if [ -s {input.sam} ]; then
            racon -t {threads} -u {input.reads} {input.sam} {input.draft} > {output.fa} 2> {log}
        else
            : > {output.fa}; echo "nothing to polish" > {log}
        fi
        """


FINAL_DRAFT = POLISH + (f"/racon{RACON_ROUNDS}.fasta" if RACON_ROUNDS > 0 else "/seeds.fasta")


if USE_MEDAKA:

    rule medaka:
        """Medaka once per sample on the Racon draft; the reads-to-draft BAM gives the read counts."""
        input:
            draft=FINAL_DRAFT,
            reads=POLISH + "/reads.fq.gz"
        output:
            fa=POLISH + "/polished.fasta",
            bam=temp(POLISH + "/reads_to_draft.bam")
        log:
            "logs/consensus/{sample}_medaka.log"
        threads: THREADS_POLISH
        params:
            outdir=POLISH + "/medaka",
            opts=(f"-m {cfg('consensus', 'medaka_model')}" if cfg("consensus", "medaka_model")
                  else ("--bacteria" if cfg("consensus", "medaka_bacteria", default=False) else "")),
            keep=int(bool(cfg("consensus", "keep_intermediate", default=False)))
        resources:
            mem_mb=mem_mb(THREADS_POLISH, 8000), runtime=RUNTIME_MIN
        container:
            CONTAINERS["medaka"]
        shell:
            r"""
            set -euo pipefail
            if [ -s {input.draft} ]; then
                rm -rf {params.outdir}
                medaka_consensus -i {input.reads} -d {input.draft} -o {params.outdir} -t {threads} {params.opts} > {log} 2>&1
                cp {params.outdir}/consensus.fasta {output.fa}
                mv {params.outdir}/calls_to_draft.bam {output.bam}
                if [ "{params.keep}" = "0" ]; then rm -f {params.outdir}/consensus_probs.hdf {params.outdir}/calls_to_draft.bam.bai; fi
            else
                : > {output.fa}; : > {output.bam}; echo "no draft" > {log}
            fi
            """

else:

    rule map_to_draft:
        """Without Medaka: the Racon draft is the consensus; map the reads once more for the read counts."""
        input:
            draft=FINAL_DRAFT,
            reads=POLISH + "/reads.fq.gz"
        output:
            fa=POLISH + "/polished.fasta",
            bam=temp(POLISH + "/reads_to_draft.bam")
        log:
            "logs/consensus/{sample}_map_to_draft.log"
        threads: THREADS_POLISH
        resources:
            mem_mb=mem_mb(THREADS_POLISH, 4000), runtime=RUNTIME_MIN
        container:
            CONTAINERS["minimap2_samtools"]
        shell:
            r"""
            set -euo pipefail
            cp {input.draft} {output.fa}
            if [ -s {input.draft} ] && [ -s {input.reads} ]; then
                minimap2 -t {threads} -ax map-ont --secondary=no {input.draft} {input.reads} 2> {log} \
                  | samtools sort -@ 2 -m 1G -o {output.bam} - 2>> {log}
            else
                : > {output.bam}; echo "no draft or no reads" > {log}
            fi
            """


rule count_mapped:
    """Primary mapped reads per consensus (reads of clusters below min_size are thereby re-assigned)."""
    input:
        bam=POLISH + "/reads_to_draft.bam"
    output:
        tsv=POLISH + "/reads_per_consensus.tsv"
    log:
        "logs/consensus/{sample}_count.log"
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["samtools"]
    shell:
        r"""
        set -euo pipefail
        if [ -s {input.bam} ]; then
            samtools view -F 0x904 {input.bam} 2> {log} | cut -f3 | sort | uniq -c | awk '{{print $2"\t"$1}}' > {output.tsv}
        else
            : > {output.tsv}; : > {log}
        fi
        """


rule consensus_collect:
    """Per amplicon: consensus.fasta and cluster_sizes.tsv (cluster size and reads mapped to each consensus)."""
    input:
        table=POLISH + "/seeds.tsv",
        polished=POLISH + "/polished.fasta",
        counts=POLISH + "/reads_per_consensus.tsv"
    output:
        cons=expand("results/consensus/{{sample}}/{amplicon}/consensus.fasta", amplicon=AMP_LIST),
        sizes=expand("results/consensus/{{sample}}/{amplicon}/cluster_sizes.tsv", amplicon=AMP_LIST)
    log:
        "logs/consensus/{sample}_collect.log"
    params:
        amplicons=AMP_LIST
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["python"]
    script:
        "../scripts/collect_consensus.py"
