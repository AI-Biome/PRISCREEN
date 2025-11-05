"""
Combine per-amplicon summaries into a sample-level summary
"""
import pandas as pd
import glob
import os

def main(snakemake):
    """Combine all amplicon summaries for a sample"""

    rows = []
    for p in snakemake.input.summaries:
        amp = p.split("/")[-2]
        if os.path.exists(p) and os.path.getsize(p) > 0:
            df = pd.read_csv(p, sep="\t")
            df.insert(0, "amplicon", amp)
            rows.append(df)

    if rows:
        all_df = pd.concat(rows, ignore_index=True)
    else:
        all_df = pd.DataFrame(
            columns=["amplicon", "consensus", "best_species", "pident", "qcov", "status"]
        )

    all_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

if __name__ == "__main__":
    main(snakemake)
