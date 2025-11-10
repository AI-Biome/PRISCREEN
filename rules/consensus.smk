"""
Rules for consensus sequence generation
"""

rule build_racon:
    output:
        generic = "software/racon_build/generic/bin/racon",
        avx2    = "software/racon_build/avx2/bin/racon"
    conda:
        "../envs/racon-build.yaml"
    threads: 4
    shell:
        r"""
        set -euo pipefail
        rm -rf software/racon_src software/racon_build
        git clone --recursive https://github.com/lbcb-sci/racon.git software/racon_src

        # Generic (portable)
        mkdir -p software/racon_build/generic
        cd software/racon_build/generic
        cmake -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_CXX_FLAGS="-O3 -march=x86-64 -mtune=generic" \
              -DCMAKE_C_FLAGS="-O3 -march=x86-64 -mtune=generic" \
              ../../racon_src
        make -j {threads}
        install -D bin/racon ../../../racon/bin/racon_generic

        # AVX2 (fast)
        cd ..
        mkdir -p avx2
        cd avx2
        cmake -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_CXX_FLAGS="-O3 -mavx2 -mtune=native" \
              -DCMAKE_C_FLAGS="-O3 -mavx2 -mtune=native" \
              ../../racon_src
        make -j {threads}
        install -D bin/racon ../../../racon/bin/racon_avx2
        """

def select_racon_output():
    outs = rules.build_racon.output
    if has_flag("avx2"):
        return outs.avx2

    return outs.generic

rule consensus_polish:
    input:
        racon = select_racon_output(),
        fq = "results/bin/{sample}/fastq/{amplicon}.fq.gz",
        cents = lambda wc: (
            f"results/cluster/{wc.sample}/{wc.amplicon}/centroids.fasta"
            if CLUSTER else
            []
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
