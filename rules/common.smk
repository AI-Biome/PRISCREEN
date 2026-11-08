"""Configuration parsing, sample sheet and helpers shared by all rule modules."""

import re


def cfg(*keys, default=None):
    """config value at a nested key path, or default when any level is missing."""
    node = config
    for k in keys:
        if not isinstance(node, dict) or k not in node or node[k] is None:
            return default
        node = node[k]
    return node


def abspath_from_config(p):
    """Paths in the config are relative to the workflow directory (where the Snakefile is)."""
    if p is None:
        return None
    return p if os.path.isabs(p) else os.path.join(workflow.basedir, p)


# ----------------------------------------------------------------------------
# Samples
# ----------------------------------------------------------------------------
SAMPLES_TSV = abspath_from_config(cfg("samples", default="config/samples.tsv"))
SAMPLES_DF = pd.read_csv(SAMPLES_TSV, sep="\t", comment="#", dtype=str).dropna(how="all")
for col in ("sample", "fastq"):
    if col not in SAMPLES_DF.columns:
        raise ValueError(f"{SAMPLES_TSV} needs the columns 'sample' and 'fastq'")
SAMPLES = SAMPLES_DF["sample"].str.strip().tolist()
if len(set(SAMPLES)) != len(SAMPLES):
    raise ValueError("duplicate sample names in " + SAMPLES_TSV)
for s in SAMPLES:
    if not s or "/" in s or " " in s:
        raise ValueError(f"invalid sample name {s!r}: no spaces or slashes")
if SAMPLES_DF["fastq"].isna().any():
    raise ValueError("empty 'fastq' cell in " + SAMPLES_TSV)
FASTQ = {r.sample.strip(): abspath_from_config(r.fastq.strip()) for r in SAMPLES_DF.itertuples()}
for s, f in FASTQ.items():
    if not os.path.exists(f):
        raise FileNotFoundError(f"FASTQ of sample {s} not found: {f}")

# ----------------------------------------------------------------------------
# Panel (PRISTINE amplicons) and rRNA amplicons
# ----------------------------------------------------------------------------
PANEL_FASTA_CFG = abspath_from_config(cfg("panel", "fasta"))
PRISTINE_DIRS = cfg("panel", "pristine_output")
if isinstance(PRISTINE_DIRS, str):
    PRISTINE_DIRS = [PRISTINE_DIRS]
PRISTINE_DIRS = [abspath_from_config(d) for d in (PRISTINE_DIRS or [])]



def check_name(name, what):
    """Amplicon names become file names, wildcards and reference-name prefixes."""
    if not name or re.search(r"[\s/|]", name):
        raise ValueError(f"invalid {what} {name!r}: must not be empty or contain spaces, '/' or '|'")


def fasta_ids(path):
    """First header token of every record."""
    with open(path) as fh:
        return [line[1:].split()[0] if len(line) > 1 and not line[1:].isspace() else "" for line in fh if line.startswith(">")]


def check_unique_ids(path, what):
    ids = fasta_ids(path)
    dup = sorted({i for i in ids if ids.count(i) > 1}) if len(set(ids)) != len(ids) else []
    if dup:
        raise ValueError(f"duplicate sequence ids in {what} {path}: {', '.join(dup[:5])}. Make every record id unique "
                         "(panel: add '|genome=<id>' to the header).")


def panel_amplicons_from_fasta(path):
    names = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].split()[0] if line[1:].strip() else ""
            if "|species=" not in head:
                raise ValueError(f"panel header must be 'amplicon|species=Genus_species[|key=value...]' without spaces: "
                                 f"{line.strip()!r} (see README: panel format, scripts/pristine_to_panel.py)")
            name = head.split("|")[0]
            check_name(name, "panel amplicon name")
            names.add(name)
    return sorted(names)


if PRISTINE_DIRS:
    # The panel is built by a rule; the amplicon list comes from the PRISTINE file names.
    PANEL_FASTA = "results/panel/panel.fasta"
    PANEL_AMPLICONS = amplicon_names(find_amplicon_files(PRISTINE_DIRS))
    for _n in PANEL_AMPLICONS:
        check_name(_n, "PRISTINE locus name")
elif PANEL_FASTA_CFG:
    if not os.path.exists(PANEL_FASTA_CFG):
        raise FileNotFoundError(f"panel.fasta not found: {PANEL_FASTA_CFG}")
    PANEL_FASTA = PANEL_FASTA_CFG
    PANEL_AMPLICONS = panel_amplicons_from_fasta(PANEL_FASTA)
    check_unique_ids(PANEL_FASTA, "panel")
else:
    PANEL_FASTA = None
    PANEL_AMPLICONS = []

RRNA = {str(k): v for k, v in (cfg("rrna", default={}) or {}).items()}
for name, spec in RRNA.items():
    if not isinstance(spec, dict) or not spec.get("database"):
        raise ValueError(f"rrna.{name} needs a 'database' FASTA")
    check_name(name, "rRNA amplicon name")
    if name == "panel":
        raise ValueError("'panel' is reserved for the PRISTINE panel database; rename the rRNA amplicon")
    if name in PANEL_AMPLICONS:
        raise ValueError(f"rRNA amplicon {name!r} clashes with a panel amplicon of the same name")
    spec["database"] = abspath_from_config(spec["database"])
    spec["bin_reference"] = abspath_from_config(spec.get("bin_reference")) or spec["database"]
    for key in ("database", "bin_reference"):
        if not os.path.exists(spec[key]):
            raise FileNotFoundError(f"rrna.{name}.{key} not found: {spec[key]}")
    check_unique_ids(spec["bin_reference"], f"rrna.{name}.bin_reference")
RRNA_NAMES = sorted(RRNA)
_bin_refs = [RRNA[n]["bin_reference"] for n in RRNA_NAMES]
if len(set(_bin_refs)) != len(_bin_refs):
    raise ValueError("two rRNA amplicons share the same bin_reference; reads could not be told apart. "
                     "Give each amplicon its own bin_reference (e.g. a bacterial and a eukaryotic subset).")

AMP_LIST = PANEL_AMPLICONS + RRNA_NAMES
if not AMP_LIST:
    raise ValueError("nothing to do: configure panel.fasta / panel.pristine_output and/or rrna amplicons")

# Databases: one for the whole panel, one per rRNA amplicon.
DATABASES = {}
if PANEL_AMPLICONS:
    DATABASES["panel"] = PANEL_FASTA
for n in RRNA_NAMES:
    DATABASES[n] = RRNA[n]["database"]


def database_of(amplicon):
    return "panel" if amplicon in PANEL_AMPLICONS else amplicon


wildcard_constraints:
    sample="|".join(map(re.escape, SAMPLES)),
    amplicon="|".join(map(re.escape, AMP_LIST)),
    db="|".join(map(re.escape, DATABASES)),

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
FILTER = bool(cfg("filter", "enable", default=False))
THREADS_MAP = int(cfg("threads", "map", default=8))
THREADS_POLISH = int(cfg("threads", "polish", default=8))
THREADS_MMSEQS = int(cfg("threads", "mmseqs", default=8))
CLUSTER_ID = float(cfg("cluster", "id_threshold", default=0.93))
CLUSTER_MIN_SIZE = int(cfg("cluster", "min_size", default=3))
IDENTIFY = bool(cfg("identify", "enable", default=True))
HAPLOTYPES = bool(cfg("haplotypes", "enable", default=True))
HAPLOTYPE_DEFAULTS = {"min_reads": 5, "min_allele_frac": 0.15, "min_depth": 20, "min_base_quality": 10,
                      "min_linkage": 0.8, "min_variants": 1, "max_haplotypes": 4,
                      "indels": True, "hp_min_run": 4, "hp_min_allele_frac": 0.2}
# Consensus sequences the identification works on: per haplotype or straight from polishing.
CONS_DIR = "results/haplotypes" if HAPLOTYPES else "results/consensus"
INDEX_RRNA_DB = bool(cfg("identify", "index_rrna_db", default=False))

MEM_MB_PER_CPU = int(cfg("resources", "mem_mb_per_cpu", default=2600))

# ----------------------------------------------------------------------------
# Containers (used with --software-deployment-method apptainer): public BioContainers,
# one per tool; overridable per key under `containers:` in the config (URI or local .sif).
# ----------------------------------------------------------------------------
CONTAINER_DEFAULTS = {
    "minimap2": "docker://quay.io/biocontainers/minimap2:2.31--h118bc1c_0",
    "minimap2_samtools": "docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:ad4c240457048f1d3cc0bd8c8fa39bf27d222eb7-0",
    "samtools": "docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1",
    "vsearch": "docker://quay.io/biocontainers/vsearch:2.32.0--hd2be7a0_0",
    "mmseqs2": "docker://quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0",
    "chopper": "docker://quay.io/biocontainers/chopper:0.15.0--h7f49ad2_0",
    "racon": "docker://quay.io/biocontainers/racon:1.5.0--h077b44d_8",
    "medaka": "docker://quay.io/biocontainers/medaka:2.2.1--py311h063c2b6_0",   # also minimap2, samtools, pysam, numpy
    "python": "docker://quay.io/biocontainers/pandas:2.2.1",                   # python + pandas for the table scripts
}
CONTAINERS = {}
for _key, _default in CONTAINER_DEFAULTS.items():
    _img = cfg("containers", _key, default=_default)
    if not str(_img).startswith(("docker://", "shub://", "library://", "oras://", "http://", "https://")) and not os.path.isabs(_img):
        _img = os.path.join(workflow.basedir, _img)
    CONTAINERS[_key] = _img
RUNTIME_MIN = int(cfg("resources", "runtime_min", default=120))


def mem_mb(threads=1, minimum=1000):
    """SLURM memory request: per-CPU share times threads (ignored in local runs)."""
    return max(minimum, MEM_MB_PER_CPU * threads)



