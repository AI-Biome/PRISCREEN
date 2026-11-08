"""Reference panel, mapping reference and MMseqs2 databases."""

import shlex


rule build_panel_from_pristine:
    """Convert PRISTINE primers/amplicons/*_amplicons.fasta into the PRISCREEN panel format."""
    input:
        PRISTINE_DIRS and find_amplicon_files(PRISTINE_DIRS)
    output:
        fa="results/panel/panel.fasta"
    log:
        "logs/panel/build_panel_from_pristine.log"
    params:
        dirs=" ".join(shlex.quote(d) for d in PRISTINE_DIRS)
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    shell:
        "python {workflow.basedir}/scripts/pristine_to_panel.py {params.dirs} -o {output.fa} 2> {log}"


rule mapping_reference:
    """One FASTA for binning: panel records as they are, rRNA bin references prefixed '<amplicon>|'."""
    input:
        panel=[PANEL_FASTA] if PANEL_AMPLICONS else [],
        rrna=[RRNA[n]["bin_reference"] for n in RRNA_NAMES]
    output:
        fa="results/panel/mapping_reference.fasta",
        fai="results/panel/mapping_reference.fasta.fai"
    log:
        "logs/mapping_reference.log"
    params:
        rrna_names=RRNA_NAMES
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    run:
        with open(output.fa, "w") as out:
            for p in input.panel:
                with open(p) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            line = ">" + line[1:].split()[0] + "\n"
                        out.write(line)
            for name, p in zip(params.rrna_names, input.rrna):
                with open(p) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            line = f">{name}|{line[1:].split()[0]}\n"
                        out.write(line)
        shell("samtools faidx {output.fa}")


rule mmseqs_createdb:
    """MMseqs2 nucleotide database for the panel or one rRNA database (indexed for large rRNA databases)."""
    input:
        fasta=lambda wc: DATABASES[wc.db]
    output:
        done="results/db/{db}/db.done"
    log:
        "logs/db/{db}.log"
    threads: THREADS_MMSEQS
    params:
        prefix="results/db/{db}/db",
        tmp="results/db/{db}/tmp",
        index=lambda wc: "1" if (wc.db != "panel" and INDEX_RRNA_DB) else "0"
    resources:
        mem_mb=mem_mb(THREADS_MMSEQS, 4000), runtime=RUNTIME_MIN
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/db/{wildcards.db}
        rm -f {params.prefix}.idx {params.prefix}.idx.index {params.prefix}.idx.dbtype
        mmseqs createdb {input.fasta} {params.prefix} --dbtype 2 -v 1 > {log} 2>&1
        if [ "{params.index}" = "1" ]; then
            mmseqs createindex {params.prefix} {params.tmp} --search-type 3 --threads {threads} -v 1 >> {log} 2>&1
            rm -rf {params.tmp}
        fi
        touch {output.done}
        """
