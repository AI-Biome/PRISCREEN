# DataProcessing

Oxford Nanopore Technology (ONT) Custom Amplicon Analysis Pipeline

Version: 0.1.0-beta.1

## Table of Contents

- [What is this?](#what-is-this)
- [Quick Start Guide](#quick-start-guide)
  - [Step 1: Install Required Software](#step-1-install-required-software)
  - [Step 2: Configure Your Analysis](#step-2-configure-your-analysis)
  - [Step 3: Run the Workflow](#step-3-run-the-workflow)
    - [Run on Local Machine](#run-on-local-machine)
    - [Run on HPC Cluster (SLURM)](#run-on-hpc-cluster-slurm)
  - [Step 4: Collect Your Results](#step-4-collect-your-results)
- [Workflow Overview](#workflow-overview)
- [Utility Scripts](#utility-scripts)
- [Advanced Configuration](#advanced-configuration)
- [Output Files](#output-files)
- [Known Issues](#known-issues)
  - [Troubleshooting Conda Environment Issues](#troubleshooting-conda-environment-issues)
  - [SLURM-Specific Troubleshooting](#slurm-specific-troubleshooting)
- [Citation](#citation)

## What is this?

This repository contains tools for processing Oxford Nanopore Technology (ONT) amplicon sequencing data from mixed samples containing both 16S rRNA and custom amplicons. The main pipeline identifies which species are present in each custom amplicon by:

1. Mapping ONT reads to a reference panel of amplicons
2. Binning reads by amplicon type
3. Clustering reads to reduce sequencing errors
4. Generating consensus sequences through polishing (Racon + Medaka)
5. Identifying species by sequence similarity search (MMseqs2)

**Use cases:**
- Analyzing environmental samples with multiple custom amplicons
- Species identification in targeted amplicon sequencing
- Quality control and validation of ONT amplicon data
- Creating synthetic test datasets with controlled error profiles

## Quick Start Guide

### Step 1: Install Required Software

You need two things:
1. **Conda/Mamba** (for managing software dependencies)
2. **Snakemake** (workflow manager)

#### Installing Snakemake (one-time setup)

```bash
# Navigate to a directory with enough space (several GB needed)
cd /path/to/your/workspace

# Download and install Miniconda (if you don't have conda/mamba)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow the prompts (accept defaults)
source ~/.bashrc  # Reload your shell

# Install mamba (faster than conda)
conda install -n base -c conda-forge mamba

# Create a Snakemake environment
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake

# Optional: Install SLURM plugin for HPC clusters
pip install snakemake-executor-plugin-slurm
```

#### Downloading This Workflow

```bash
git clone https://github.com/YOUR_USERNAME/DataProcessing.git
cd DataProcessing
```

**Note**: All bioinformatics tools (minimap2, samtools, MMseqs2, vsearch, racon, medaka, seqtk) are automatically installed via conda environments - you don't need to install them manually.

### Step 2: Configure Your Analysis

#### Prepare Reference Panel

Create a FASTA file with your reference amplicons in the format:

```
>amplicon_id|species=species_name
SEQUENCE
```

Example (data/reference_amplicons.fasta):

```
>amp1|species=SpA
ACGTACGTACGTACGT...
>amp1|species=SpB
ACGTACGTTTGTACGT...
>amp2|species=SpA
GCTAGCTAGCTAGCTA...
>amp2|species=SpC
GCTAGCTAGCCCGCTA...
```

**Important**: Each amplicon can have multiple species variants. The pipeline will:
1. Map reads to all variants
2. Extract reads for each amplicon type (amp1, amp2, etc.)
3. Generate consensus sequences
4. Identify which species match best

#### Configure Samples

Edit config/samples.tsv to list your samples:

```tsv
sample	fastq
S1	data/S1.fastq.gz
S2	data/S2.fastq.gz
```

**Columns:**
- sample: Sample name (no spaces)
- fastq: Path to FASTQ file (can be gzipped, relative to Snakefile location)

#### Adjust Pipeline Parameters

Edit config/config.yaml to customize the pipeline. See [Advanced Configuration](#advanced-configuration) for details.

### Step 3: Run the Workflow

#### Test Your Configuration

Before running, check that everything is configured correctly:

```bash
# Activate Snakemake environment
mamba activate snakemake

# Dry run - shows what would be executed
snakemake --cores 1 --use-conda --dry-run
```

#### Run on Local Machine

```bash
# Activate Snakemake environment
mamba activate snakemake

# Run pipeline using 8 cores
snakemake --cores 8 --use-conda --conda-frontend conda
```

Replace 8 with the number of CPU cores available on your machine.

**Note:** The `--conda-frontend conda` flag is recommended if you encounter mamba/conda compatibility issues. It forces Snakemake to use conda instead of mamba for environment creation, which can prevent ImportError issues with numpy/pandas.

#### Run on HPC Cluster (SLURM)

The pipeline includes SLURM configuration for running on HPC clusters. Resource allocation (CPUs, memory, runtime) is configured via `config.ini` and automatically applied to all rules.

**Method 1: Using SLURM Profile (Recommended)**

```bash
# Activate Snakemake environment
mamba activate snakemake

# Submit jobs to SLURM using the pre-configured profile
snakemake --profile profiles/slurm
```

**Method 2: Command-line Options**

```bash
# Submit jobs to SLURM with explicit options
snakemake \
  --executor slurm \
  --jobs 20 \
  --use-conda \
  --default-resources slurm_partition=batch slurm_account=none
```

**Command explanation:**
- `--executor slurm`: Use SLURM job scheduler
- `--jobs 20`: Run up to 20 jobs simultaneously
- `--use-conda`: Use conda environments for tools
- `--default-resources`: Set default SLURM partition and account

**Configuring SLURM Resources**

Edit `config.ini` to adjust resource allocation for your cluster:

```ini
[SLURM_ARGS]
cpus_per_task = 48          # CPUs per SLURM task
mem_of_node = 126000        # Total node memory in MB
max_runtime = 4320          # Max runtime in minutes (72h)
slurm_partition = batch     # SLURM partition name
```

**Parameter explanation:**
- `cpus_per_task`: Number of CPUs to request per job
- `mem_of_node`: Total memory available on the node in MB (memory per job = mem_of_node / cpus_per_task)
- `max_runtime`: Maximum runtime in minutes
- `slurm_partition`: Name of the SLURM partition to use

Edit `profiles/slurm/config.yaml` to adjust cluster-specific settings:

```yaml
executor: slurm
jobs: 20                    # Max simultaneous jobs
default-resources:
  - slurm_account=none      # Your SLURM account
  - slurm_partition=batch   # Your partition name
```

#### Monitor Progress

```bash
# Check SLURM queue (if using cluster)
squeue -u $USER

# Watch jobs in real-time
watch -n 5 'squeue -u $USER'

# Check for output files
ls -lh results/identify/*/species_summary.tsv

# View Snakemake log
tail -f .snakemake/log/*.log

# View SLURM logs for specific rules
tail -f .snakemake/slurm_logs/rule_*/slurm-*.out
```

### Step 4: Collect Your Results

Once the workflow completes, find your results in results/identify/{sample}/:

**Main output file:**

```
results/identify/{sample}/species_summary.tsv
```

This file contains species identifications for all amplicons in the sample:

```tsv
amplicon	consensus	best_species	pident	qcov	status
amp1	consensus_1	SpA	0.995	0.98	OK
amp1	consensus_2	SpB	0.987	0.96	OK
amp2	consensus_1	SpC	0.992	0.99	AMBIGUOUS
```

**Columns:**
- amplicon: Amplicon identifier
- consensus: Consensus sequence identifier
- best_species: Top matching species
- pident: Percent identity (0-1 scale)
- qcov: Query coverage (0-1 scale)
- status: OK, AMBIGUOUS, or NO_HIT

**Status meanings:**
- OK: Clear species identification
- AMBIGUOUS: Multiple species match within top_delta threshold
- NO_HIT: No species match above quality thresholds

## Workflow Overview

This Snakemake pipeline processes ONT amplicon data to identify species in mixed samples. The pipeline uses conda environments to manage all dependencies automatically.

### Pipeline Steps

1. mmseqs_createdb - Build searchable database from reference panel
2. map_to_panel - Align reads to reference amplicons (minimap2 -ax map-ont)
3. amplicon_readnames - Extract read names for each amplicon type
4. amplicon_fastq - Extract FASTQ records for individual amplicons (seqtk)
5. cluster_vsearch - Cluster reads by similarity (vsearch, default 99.5% identity)
6. consensus_polish - Generate consensus sequences:
   - Optional: Use cluster centroids as seeds (if clustering enabled)
   - Fallback: Use first read as seed (if clustering disabled)
   - Racon polishing (default: 2 rounds)
   - Medaka error correction
7. mmseqs_search - Search consensus sequences against reference database
8. summarize_hits - Identify best species match with quality filtering
9. sample_summary - Aggregate results per sample across all amplicons

**Total runtime**: Varies by sample size and number of amplicons
- Small samples (< 10K reads): ~30 minutes
- Medium samples (10K-100K reads): 1-3 hours
- Large samples (> 100K reads): 3-6 hours

## Utility Scripts

The repository includes standalone utility scripts in the scripts/ directory:

### FASTQ Error Simulation

**Purpose**: Introduce controlled sequencing errors into FASTQ files for testing.

**Usage:**

```bash
python scripts/fastq_error_simulation.py \
  --input data/input.fastq.gz \
  --output data/output.fastq.gz \
  --error_rate 0.01 \
  --seed 42
```

**Parameters:**
- --input: Input FASTQ file (plain or gzipped)
- --output: Output FASTQ file (same format as input)
- --error_rate: Per-base substitution error rate (0-1), e.g., 0.01 = 1%
- --seed: Random seed for reproducibility (optional)
- --low_quality_char: Quality score for mutated bases (default: "!")
- --preserve_quality: Keep original quality scores even for mutated bases

**Features:**
- Case-preserving mutations (maintains upper/lowercase)
- Supports gzipped files
- Optionally updates quality scores

### FASTQ Random Sampling

**Purpose**: Subsample a specified number of reads from a FASTQ file.

**Usage:**

```bash
python scripts/fastq_random_sample.py \
  --input data/input.fastq \
  --output data/output.fastq \
  --num_reads 10000 \
  --seed 42
```

**Parameters:**
- --input: Input FASTQ file (plain text only, no gzip support)
- --output: Output FASTQ file
- --num_reads: Number of reads to sample
- --seed: Random seed for reproducibility (optional)

**Note**: This script loads the entire FASTQ into memory. Not recommended for very large files (> 1GB).

### Read Mapping to Reference (INCOMPLETE)

**Warning**: This script is incomplete and not functional.

**File**: scripts/map_fastq_to_fasta.py

**Intended purpose**: Map query reads to target sequences using minimap2 or pblat.

**Status**: Contains multiple bugs and incomplete implementations. See ISSUES.md for details.

## Advanced Configuration

Edit config/config.yaml to customize the pipeline:

### Clustering Parameters

```yaml
cluster:
  enable: true              # Set to false to skip clustering
  id_threshold: 0.995       # Clustering identity (0.99-1.0 typical)
```

**When to disable clustering:**
- Very low coverage (< 10 reads per amplicon)
- High-accuracy basecalling already performed
- Testing/debugging

### Consensus Polishing

```yaml
consensus:
  racon_rounds: 2           # More rounds = better accuracy but slower
  medaka_model: "r1041_e82_400bps_hac_v4.2.0"
```

**Medaka models** (choose based on your sequencing chemistry):
- r1041_e82_400bps_hac_v4.2.0 - R10.4.1, 400 bps, HAC basecalling
- r1041_e82_400bps_sup_v4.2.0 - R10.4.1, 400 bps, SUP basecalling
- r941_min_hac_g507 - R9.4.1, MinION, HAC basecalling

See Medaka documentation for full list: https://github.com/nanoporetech/medaka#models

### Species Identification Thresholds

```yaml
identify:
  min_pident: 0.985         # Minimum identity (0-1 scale)
  min_qcov: 0.90            # Minimum query coverage (0-1 scale)
  top_delta: 0.01           # Ambiguity threshold
```

**Threshold guidelines:**
- min_pident: 0.97-0.99 for closely related species, 0.90-0.95 for genus-level
- min_qcov: 0.80-0.95 typical, depends on amplicon length variation
- top_delta: 0.01-0.02 for strict matches, 0.05 for more permissive

**Status determination:**
- OK: Top hit passes thresholds, second-best hit differs by > top_delta
- AMBIGUOUS: Multiple hits within top_delta of top hit
- NO_HIT: No hits pass min_pident and min_qcov

### Thread Configuration

```yaml
threads:
  map: 8                    # Minimap2 alignment threads
  polish: 8                 # Racon and Medaka threads
  mmseqs: 8                 # MMseqs2 search threads
```

Adjust based on available CPU cores. More threads = faster processing but higher memory usage.

## Output Files

### Per-Sample Summary

**File**: results/identify/{sample}/species_summary.tsv

Tab-delimited file with species identifications for all amplicons.

### Per-Amplicon Detailed Results

**File**: results/identify/{sample}/{amplicon}/summary.tsv

Detailed results for each amplicon separately.

### Raw MMseqs2 Output

**File**: results/identify/{sample}/{amplicon}/mmseqs_hits.tsv

Raw similarity search results (13-column format).

### Consensus Sequences

**File**: results/consensus/{sample}/{amplicon}/consensus.fasta

Polished consensus sequences after Racon + Medaka.

### Read Bins

**File**: results/bin/{sample}/fastq/{amplicon}.fq.gz

Reads binned by amplicon type.

## Known Issues

See ISSUES.md for detailed documentation of known bugs and limitations.

**Critical issues:**
- scripts/map_fastq_to_fasta.py is non-functional (missing imports, incomplete logic)
- Snakefile line 246 has incorrect percentage comparison logic

**Design inconsistencies:**
- Inconsistent file format support across utility scripts
- Memory inefficiency in scripts/fastq_random_sample.py

### Troubleshooting Conda Environment Issues

If you encounter conda dependency conflicts when creating environments (especially with `medaka`, `samtools`, or `htslib`), try setting channel priority to flexible:

```bash
conda config --set channel_priority flexible
```

This allows conda to search across all configured channels to find compatible package versions, rather than strictly prioritizing packages from specific channels. This is particularly important for bioinformatics workflows that mix packages from `bioconda` and `conda-forge`.

### SLURM-Specific Troubleshooting

**Jobs failing with memory errors:**
- Increase `mem_of_node` in `config.ini`
- Or reduce `cpus_per_task` (allocates more memory per job: mem_of_node / cpus_per_task)

**Jobs timing out:**
- Increase `max_runtime` in `config.ini`
- Check partition time limits: `sinfo -o "%P %l"`

**Permission denied on partition:**
- Check available partitions: `sinfo`
- Update `slurm_partition` in `profiles/slurm/config.yaml`
- Verify access: `sacctmgr show user $USER withassoc`

**Too many jobs queued:**
- Reduce `jobs` value in `profiles/slurm/config.yaml`
- Check cluster fair-use policies

**Conda environments not building on compute nodes:**
- Pre-create environments: `snakemake --use-conda --conda-create-envs-only --cores 1`
- Check internet connectivity on compute nodes
- Consider using `--conda-frontend conda` instead of mamba

## Citation

If you use this pipeline, please cite the tools it depends on:

- Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100.
- Samtools: Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008.
- Racon: Vaser, R., et al. (2017). Fast and accurate de novo genome assembly from long uncorrected reads. Genome Research, 27(5), 737-746.
- Medaka: Oxford Nanopore Technologies. Medaka: Sequence correction provided by ONT Research. https://github.com/nanoporetech/medaka
- MMseqs2: Steinegger, M., & Soding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35(11), 1026-1028.
- vsearch: Rognes, T., et al. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.
- Snakemake: Molder, F., et al. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, 33
