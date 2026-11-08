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
        dirs=PRISTINE_DIRS
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["python"]
    script:
        "../scripts/build_panel.py"


rule mapping_reference:
    """One FASTA for binning: panel records as they are (id = first header token), rRNA bin references
    prefixed '<amplicon>|'; indexed with samtools faidx."""
    input:
        panel=[PANEL_FASTA] if PANEL_AMPLICONS else [],
        rrna=[RRNA[n]["bin_reference"] for n in RRNA_NAMES]
    output:
        fa="results/panel/mapping_reference.fasta",
        fai="results/panel/mapping_reference.fasta.fai"
    log:
        "logs/mapping_reference.log"
    params:
        panel=lambda wc, input: " ".join(shlex.quote(p) for p in input.panel),
        rrna=lambda wc, input: " ".join(shlex.quote(f"{n}={p}") for n, p in zip(RRNA_NAMES, input.rrna))
    resources:
        mem_mb=mem_mb(1), runtime=RUNTIME_MIN
    container:
        CONTAINERS["samtools"]
    shell:
        r"""
        set -euo pipefail
        : > {output.fa}
        for f in {params.panel} ""; do
            [ -n "$f" ] || continue
            awk '/^>/ {{print ">" substr($1, 2); next}} {{print}}' "$f" >> {output.fa}
        done
        for np in {params.rrna} ""; do
            [ -n "$np" ] || continue
            n="${{np%%=*}}"; p="${{np#*=}}"
            awk -v N="$n" '/^>/ {{print ">" N "|" substr($1, 2); next}} {{print}}' "$p" >> {output.fa}
        done
        samtools faidx {output.fa} 2> {log}
        """


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
    container:
        CONTAINERS["mmseqs2"]
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
