# Known Issues and Bugs

This document details known bugs, incomplete implementations, and design inconsistencies in the DataProcessing repository.

Last updated: 2025-11-05

## Critical Issues

### FIXED Issue 1: map_fastq_to_fasta.py is Non-Functional

**File**: scripts/map_fastq_to_fasta.py

**Severity**: Critical - Script cannot run

**Description**: The script contains multiple bugs that prevent it from executing:

#### FIXED Missing Import Statements

**Location**: Lines 8-12

**Problem**: The `re` module is used but not imported (possibly Katharina already added the import during reading).

```python
def extract_species(header):
    pattern = r">\S+\s+([A-Z][a-z]+\s[a-z]+)"
    match = re.search(pattern, header)  # ERROR: 're' is not imported
    if match:
        return match.group(1)
    return None
```

**Fix Required**: Add `import re` at the top of the file.

---

**Location**: Line 60

**Problem**: The `defaultdict` class is used but not imported from collections.

```python
query_to_targets = defaultdict(set)  # ERROR: defaultdict not imported
```

**Fix Required**: Add `from collections import defaultdict` at the top of the file.

---

#### FIXED Incomplete Implementation

**Location**: Lines 63-64

**Problem**: The minimap2 processing logic is not implemented, but the script accepts minimap2 as a valid option.

```python
if software == "minimap2":
    pass # todo  # ERROR: Not implemented
```

**Impact**: Users who select minimap2 will get no output, and the `total_mappers` count will always be 0.

**Fix Required**: Implement SAM/BAM parsing for minimap2 output, similar to the pblat implementation.

---

#### FIXED Logic Error: total_mappers Never Incremented

**Location**: Lines 61, 87

**Problem**: The `total_mappers` variable is initialized to 0 but never updated during processing.

```python
total_mappers = 0
# ... processing happens ...
csv_writer.writerow([query_file, target_file, unique_mappers, total_mappers])
# BUG: total_mappers is always 0
```

**Expected behavior**: Should be set to `len(query_to_targets)` to reflect the number of queries with at least one mapping.

**Fix Required**:
```python
total_mappers = len(query_to_targets)
csv_writer.writerow([query_file, target_file, unique_mappers, total_mappers])
```

---

#### FIXED Missing Function Arguments

**Location**: Lines 14, 99-107

**Problem**: Function signature expects `min_identity` and `max_size_diff` parameters, but the main() function doesn't pass them.

**Function definition (line 14)**:
```python
def map_reads(target_folder, query_folder, threads, software, output_dir, min_identity, max_size_diff):
```

**Function call (lines 101-107)**:
```python
summary_file = map_reads(
    target_folder=args.target_folder,
    query_folder=args.query_folder,
    threads=args.threads,
    software=args.software,
    output_dir=args.output_dir
    # ERROR: Missing min_identity and max_size_diff
)
```

**Impact**: Script will crash with TypeError when run.

**Fix Required**: Either:
1. Add `min_identity` and `max_size_diff` arguments to argparse and pass them, OR
2. Remove these parameters from the function signature and hardcode default values

---

### FIXED Issue 2: Snakefile Mathematical Bug in Ambiguity Detection

**File**: Snakefile (root directory)

**Severity**: High - Incorrect results

**Location**: Line 246

**Problem**: Incorrect percentage comparison logic when detecting ambiguous species matches.

```python
if p1 > 0 and (p1 - p2)/p1 < params.top_delta*100:
    status = "AMBIGUOUS"
```

**Context**:
- `params.top_delta` is a fraction (0.01 from config.yaml)
- `p1` and `p2` are percentages (0-100 scale from MMseqs2)
- The multiplication by 100 is incorrect

**Example of bug**:
- Config setting: `top_delta: 0.01` (1% difference allowed)
- Top hit: `p1 = 99.0`
- Second hit: `p2 = 98.5`
- Current calculation: `(99.0 - 98.5) / 99.0 = 0.00505 < 1.0` → AMBIGUOUS (wrong)
- Correct calculation: `(99.0 - 98.5) / 99.0 = 0.00505 < 0.01` → NOT ambiguous

**Impact**: Nearly all hits with minor differences will be incorrectly flagged as "AMBIGUOUS".

**Fix Required**: Remove the `*100` from the comparison:

```python
if p1 > 0 and (p1 - p2)/p1 < params.top_delta:
    status = "AMBIGUOUS"
```

**Related code**: Line 234 correctly converts config fractions to percentages for filtering:
```python
df_f = df[(df["pident"] >= params.min_pid*100) & (df["qcov"] >= params.min_qcov*100)].copy()
```

This shows the inconsistency - filtering uses `*100` correctly, but ambiguity detection does not.

---

## Design Inconsistencies

### FIXED Issue 3: Inconsistent File Format Support

**Files Affected**:
- scripts/fastq_error_simulation.py
- scripts/fastq_random_sample.py
- scripts/map_fastq_to_fasta.py

**Severity**: Medium - User experience issue

**Problem**: The three utility scripts have inconsistent support for gzipped files:

| Script | Plain FASTQ | Gzipped FASTQ |
|--------|-------------|---------------|
| fastq_error_simulation.py | Yes | Yes |
| fastq_random_sample.py | Yes | Yes |
| map_fastq_to_fasta.py | Yes (assumed) | Yes |

**Impact**: Users must decompress files before using certain scripts, creating unnecessary workflow steps.

**Recommendation**: Standardize all scripts to support both plain and gzipped FASTQ using the pattern from scripts/fastq_error_simulation.py:

```python
def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)
```

---

### FIXED Issue 4: Memory Inefficiency in fastq_random_sample.py

**File**: scripts/fastq_random_sample.py

**Severity**: Medium - Performance issue

**Location**: Lines 18-19

**Problem**: The script loads the entire FASTQ file into memory.

```python
with open(input_fastq, 'r') as f:
    lines = f.readlines()  # Loads entire file into RAM
```

**Impact**:
- Cannot process large ONT files (often 10GB+)
- Will crash with MemoryError on files > available RAM
- Inconsistent with the streaming approach used in fastq_error_simulation.py

**Comparison**:
- scripts/fastq_error_simulation.py: Streams data, can handle unlimited file size
- scripts/fastq_random_sample.py: Loads all data, limited by RAM

**Recommendation**: Implement reservoir sampling algorithm for memory-efficient random sampling:

1. Use two-pass approach: count records first, then sample
2. Alternatively, use reservoir sampling for single-pass unknown-size sampling
3. Add gzip support while refactoring

---

## Minor Issues

### FIXED Issue 5: File Extension Limitations in map_fastq_to_fasta.py

**File**: scripts/map_fastq_to_fasta.py

**Severity**: Low - Usability issue

**Location**: Lines 27, 32-33

**Problem**: Script only recognizes specific file extensions.

```python
if not target_file.endswith(".fasta"):
    continue
# ...
if not query_file.endswith(".fastq"):
    continue
```

**Impact**: Files with other common extensions are ignored:
- .fa, .fna (FASTA variants)
- .fq, .fastq.gz, .fq.gz (FASTQ variants)

**Recommendation**: Accept multiple extensions:
```python
if not target_file.endswith((".fasta", ".fa", ".fna")):
    continue
if not query_file.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
    continue
```

---

### Issue 6: Hardcoded File Extensions in Snakemake Workflow

**File**: Snakefile (root directory)

**Severity**: Low - Configuration limitation

**Location**: Throughout (wildcards and rules)

**Problem**: File paths and extensions are hardcoded, limiting flexibility.

**Examples**:
- BAM files are always named `{sample}.panel.sorted.bam`
- Output always goes to `results/` directory
- Cannot easily change intermediate file locations

**Recommendation**: Add configuration options for:
- Output directory prefix
- Temporary file directory
- Custom file naming schemes

---

## Documentation Gaps

### Issue 7: Missing Error Handling Documentation

**Severity**: Low - Documentation issue

**Problem**: No documentation exists for common error scenarios and their solutions.

**Missing topics**:
- What to do when Medaka model is not found
- How to handle empty amplicon bins (no reads map to amplicon)
- Debugging consensus polishing failures
- Handling ambiguous species identifications

**Recommendation**: Add a TROUBLESHOOTING.md file with:
- Common error messages and solutions
- How to check log files
- Debugging strategies
- FAQ section

---

### Issue 8: No Input Validation

**Files Affected**: All scripts

**Severity**: Low - Usability issue

**Problem**: Scripts do not validate inputs before processing:
- No check if reference FASTA uses correct header format
- No validation that FASTQ files are well-formed
- No check if amplicon IDs in reference match expected format

**Impact**: Cryptic errors occur deep in the pipeline when inputs are malformed.

**Recommendation**: Add validation functions:
1. Check reference FASTA headers match pattern: `>amplicon_id|species=species_name`
2. Validate FASTQ format before processing
3. Verify sample file paths exist before starting pipeline
4. Add `--validate` flag to check inputs without running pipeline

---

## Non-Issues (Intentional Design Choices)

### nanofilt Configuration Present but Disabled

**File**: ONT_custom_amplicon_analysis/config/config.yaml

**Lines**: 4-8

**Status**: This is intentional, not a bug.

```yaml
nanofilt:
  enable: false
  min_len: 200
  max_len: 0
  min_q: 9
```

The configuration exists for future use but is not currently implemented in the Snakefile. This is acceptable for development.

---

## Priority Summary

**DONE Immediate action required (Critical)**:
1. Fix scripts/map_fastq_to_fasta.py imports and implementation
2. Fix Snakefile line 246 ambiguity detection logic

**DONE Should fix soon (High priority)**:
3. Standardize gzip support across scripts
4. Optimize scripts/fastq_random_sample.py memory usage

**Nice to have (Medium priority)**:
5. DONE Improve file extension handling
6. Add input validation
7. Create troubleshooting documentation

**Low priority (Future enhancement)**:
8. Make Snakefile paths configurable
9. Implement nanofilt integration if needed

---

## Testing Recommendations

After fixing critical issues, test with:

1. Small test dataset (< 1000 reads)
2. Edge cases:
   - Empty amplicon bins
   - Single read per amplicon
   - 100% identical species sequences
   - Species with identity exactly at threshold
3. Performance test with large file (> 1GB FASTQ)
4. Verify ambiguity detection with controlled test cases

---

## Version Tracking

This issues document corresponds to the state of the repository at commit: b6331e5

Future versions should update this document as issues are resolved.
