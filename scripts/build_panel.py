#!/usr/bin/env python3
"""Snakemake wrapper around pristine_to_panel.py: build results/panel/panel.fasta from PRISTINE output."""

import sys


def main(smk):
    sys.path.insert(0, smk.scriptdir)
    from pristine_to_panel import find_amplicon_files, panel_records, write_panel

    files = find_amplicon_files(list(smk.params.dirs))
    records = panel_records(files)
    with open(smk.output.fa, "w") as out:
        n = write_panel(records, out)
    with open(smk.log[0], "w") as log:
        loci = {h.split("|")[0] for h, _ in records}
        log.write(f"{n} panel sequences for {len(loci)} amplicon(s) from {len(files)} PRISTINE file(s)\n")


if __name__ == "__main__":
    main(snakemake)  # noqa: F821 - injected by Snakemake
