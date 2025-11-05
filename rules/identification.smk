"""
Species identification and summarization rules
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

        mmseqs easy-search {input.cons} {input.db} {output.hits} tmp_mmseqs_{wildcards.sample}_{wildcards.amplicon} --threads {threads}
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
    run:
        import pandas as pd
        import os
        cols = ["qseqid","sseqid","pident","qcov","alnlen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bits"]
        df = pd.read_csv(input.hits, sep="\t", header=None, names=cols)
        if df.empty:
            out = pd.DataFrame(columns=["consensus","best_species","pident","qcov","status"])
            out.to_csv(output.tsv, sep="\t", index=False)
            return

        def parse_species(s):
            return s.split("|species=")[-1] if "|species=" in s else "UNK"
        df["species"] = df["sseqid"].map(parse_species)


        df_f = df[(df["pident"] >= params.min_pid*100) & (df["qcov"] >= params.min_qcov*100)].copy()

        out_rows = []
        for q, sub in df_f.groupby("qseqid"):
            sub = sub.sort_values(["pident","qcov","bits"], ascending=[False, False, False]).reset_index(drop=True)
            if len(sub)==0:
                out_rows.append((q, "NO_HIT", None, None, "NO_HIT"))
                continue
            top = sub.iloc[0]
            status = "OK"
            if len(sub) > 1:
                p1 = top["pident"]; p2 = sub.iloc[1]["pident"]
                if p1 > 0 and (p1 - p2)/p1 < params.top_delta:
                    status = "AMBIGUOUS"
            out_rows.append((q, top["species"], float(top["pident"])/100.0, float(top["qcov"])/100.0, status))

        out = pd.DataFrame(out_rows, columns=["consensus","best_species","pident","qcov","status"])
        out.to_csv(output.tsv, sep="\t", index=False)

rule sample_summary:
    input:
        summaries = expand("results/identify/{{sample}}/{amplicon}/summary.tsv", amplicon=AMP_LIST)
    output:
        tsv = "results/identify/{sample}/species_summary.tsv"
    conda:
        "../envs/samtools.yaml"
    run:
        import pandas as pd, glob, os
        rows = []
        for p in input.summaries:
            amp = p.split("/")[-2]
            if os.path.exists(p) and os.path.getsize(p)>0:
                df = pd.read_csv(p, sep="\t")
                df.insert(0,"amplicon",amp)
                rows.append(df)
        if rows:
            all_df = pd.concat(rows, ignore_index=True)
        else:
            all_df = pd.DataFrame(columns=["amplicon","consensus","best_species","pident","qcov","status"])
        all_df.to_csv(output.tsv, sep="\t", index=False)
