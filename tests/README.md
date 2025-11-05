# DataProcessing Test Suite

This directory contains unit tests for the DataProcessing repository scripts.

## Test Structure

```
tests/
├── README.md                          # This file
├── run_tests.sh                       # Test runner script
├── test_fastq_random_sample.py        # Tests for fastq_random_sample.py
├── test_fastq_error_simulation.py     # Tests for fastq_error_simulation.py
├── test_map_fastq_to_fasta.py         # Tests for map_fastq_to_fasta.py
└── test_snakemake_workflow.py         # Integration tests for Snakemake workflow
```

## Requirements

The test suite uses the existing `test_data/` directory in the project root:
- `test_data/data/` - Test FASTQ files (Test1.fastq.gz, Test2.fastq.gz)
- `test_data/resources/panel/` - Test reference panel (amplicon_panel.fasta)

## Running Tests

### Option 1: Using the test runner script (recommended)

```bash
cd DataProcessing
./tests/run_tests.sh
```

### Option 2: Using pytest (if available)

```bash
cd DataProcessing
python -m pytest tests/ -v
```

### Option 3: Using unittest directly

```bash
cd DataProcessing
python -m unittest discover -s tests -p "test_*.py" -v
```

### Option 4: Run individual test files

```bash
cd DataProcessing
python tests/test_fastq_random_sample.py
python tests/test_fastq_error_simulation.py
python tests/test_map_fastq_to_fasta.py
```

## Test Coverage

### test_fastq_random_sample.py

Tests for the FASTQ random sampling utility:
- ✅ Basic subsampling functionality
- ✅ Reproducibility with seeds
- ✅ Different seeds produce different results
- ✅ Plain (non-gzipped) FASTQ support
- ✅ Error handling for invalid inputs
- ✅ Gzip file detection

**9 tests total**

### test_fastq_error_simulation.py

Tests for the FASTQ error simulation utility:
- ✅ Base mutation functionality
- ✅ Case preservation in mutations
- ✅ Zero error rate produces identical output
- ✅ Reproducibility with seeds
- ✅ Mutations are introduced correctly
- ✅ Quality score updates
- ✅ Quality preservation mode
- ✅ Plain FASTQ support
- ✅ Invalid input validation
- ✅ Malformed FASTQ detection

**11 tests total**

### test_map_fastq_to_fasta.py

Tests for the read mapping utility:
- ✅ File extension detection (FASTQ/FASTA)
- ✅ Filename stem extraction
- ✅ Species name extraction from headers
- ✅ CIGAR string parsing (all operations)
- ✅ SAM optional tag parsing
- ✅ SAM flag interpretation (unmapped, secondary, supplementary)
- ✅ Integration test with minimap2 (if available)

**18 tests total**

### test_snakemake_workflow.py

Integration tests for the complete Snakemake workflow:
- ✅ Snakefile syntax validation
- ✅ Workflow DAG generation
- ✅ Rule listing and verification
- ✅ Dry-run of complete workflow
- ✅ Individual rule testing (mmseqs_createdb, map_to_panel)
- ✅ Complete workflow execution (single sample)
- ✅ Output file structure validation
- ✅ Workflow idempotence (no unnecessary re-runs)
- ✅ Modular rule file accessibility
- ✅ Configuration file validation
- ✅ Samples file validation

**13+ tests total**

## Test Examples

### Quick test run
```bash
# Run all tests
./tests/run_tests.sh
```

### Verbose output
```bash
# Using unittest
python -m unittest discover -s tests -p "test_*.py" -v

# Using pytest
python -m pytest tests/ -v
```

### Run specific test
```bash
# Test just the random sampling
python -m unittest tests.test_fastq_random_sample.TestFastqRandomSample.test_subsample_basic

# Test just the error simulation
python -m unittest tests.test_fastq_error_simulation.TestFastqErrorSimulation.test_simulate_errors_zero_rate
```

## Integration Tests

Some tests require external tools:
- **minimap2**: Required for `test_map_reads_minimap2()` integration test
  - Test will be skipped if minimap2 is not available
  - To run: ensure minimap2 is in PATH or activate appropriate conda environment

### Running Snakemake Workflow Tests

The Snakemake workflow tests require additional dependencies:

```bash
# Install snakemake (if not already installed)
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake

# Activate environment
conda activate snakemake

# Run only Snakemake tests
python -m pytest tests/test_snakemake_workflow.py -v

# Or run directly
python tests/test_snakemake_workflow.py
```

**Note**: Snakemake workflow tests can take longer to run as they execute the full pipeline. They use test data from `test_data/` with reduced parameters for faster execution.

## Continuous Integration

To run tests in CI/CD pipeline:

```bash
#!/bin/bash
# Install dependencies
conda env create -f envs/minimap2.yaml

# Run tests
conda activate minimap2
python -m pytest tests/ -v --tb=short

# Or without conda:
python -m unittest discover -s tests -p "test_*.py"
```

## Writing New Tests

To add new tests:

1. Create a new file: `tests/test_<module_name>.py`
2. Import unittest: `import unittest`
3. Add the scripts directory to path:
   ```python
   import sys
   import os
   sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
   ```
4. Create test class inheriting from `unittest.TestCase`
5. Write test methods starting with `test_`
6. Run with `./tests/run_tests.sh`

Example:
```python
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from my_script import my_function

class TestMyScript(unittest.TestCase):
    def test_my_function(self):
        result = my_function(input_value)
        self.assertEqual(result, expected_value)

if __name__ == '__main__':
    unittest.main()
```

## Test Data

The test suite uses real bioinformatics data from `test_data/`:
- **Test1.fastq.gz**: ONT reads with amplicon identifiers
- **Test2.fastq.gz**: Additional test reads
- **amplicon_panel.fasta**: Reference sequences for two species (Species_A, Species_B) across two amplicons (amp1, amp2)

The test data contains:
- Realistic ONT read quality scores
- Amplicon-based read naming convention
- Multi-species reference panel
- Gzipped FASTQ format

## Troubleshooting

### Tests fail with "test_data not found"
Ensure you're running tests from the project root:
```bash
cd DataProcessing
./tests/run_tests.sh
```

### Integration tests are skipped
Install required tools (e.g., minimap2):
```bash
conda env create -f envs/minimap2.yaml
conda activate minimap2
./tests/run_tests.sh
```

### Import errors
Ensure scripts are in the correct location:
```bash
ls -l scripts/fastq_random_sample.py
ls -l scripts/fastq_error_simulation.py
ls -l scripts/map_fastq_to_fasta.py
```

## Snakemake Workflow Testing

The workflow tests verify:

1. **Syntax and Structure**
   - Valid Snakefile syntax
   - DAG generation
   - Rule accessibility from included files

2. **Individual Rules**
   - Database creation (mmseqs_createdb)
   - Read mapping (map_to_panel)
   - Each rule can run independently

3. **Complete Workflow**
   - End-to-end execution for single sample
   - Output file generation and structure
   - Workflow idempotence (reruns don't redo work)

4. **Configuration**
   - Valid YAML config files
   - Proper samples.tsv format
   - Referenced files exist

### Test Configuration

Tests use a dedicated config at `test_data/config/config.yaml` with:
- Reduced thread counts (2 instead of 8)
- Fewer racon rounds (1 instead of 2)
- Relaxed thresholds for test data
- References to test_data paths

## Future Test Additions

Potential areas for additional tests:
- Performance benchmarks for large files
- Edge cases for consensus polishing
- MMseqs2 search validation
- Species identification accuracy
- Multi-sample workflow testing
