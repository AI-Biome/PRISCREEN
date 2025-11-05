# Known Issues and Bugs

This document details known bugs, incomplete implementations, and design inconsistencies in the DataProcessing repository.

Last updated: 2025-11-05

## Minor Issues

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
