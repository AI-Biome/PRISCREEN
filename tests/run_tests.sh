#!/bin/bash
# Test runner script for DataProcessing repository

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "========================================="
echo "DataProcessing Test Suite"
echo "========================================="
echo ""

# Check if test data exists
if [ ! -d "$PROJECT_ROOT/test_data" ]; then
    echo "ERROR: test_data directory not found at $PROJECT_ROOT/test_data"
    echo "Please ensure test data is available before running tests."
    exit 1
fi

# Activate conda environment if specified
if [ -n "$CONDA_ENV" ]; then
    echo "Activating conda environment: $CONDA_ENV"
    conda activate "$CONDA_ENV"
fi

# Run unit tests
echo "Running unit tests..."
echo ""

cd "$PROJECT_ROOT"

# Run all tests with verbose output
python -m pytest tests/ -v --tb=short 2>/dev/null || {
    # Fallback to unittest if pytest not available
    echo "pytest not found, using unittest..."
    echo ""

    python -m unittest discover -s tests -p "test_*.py" -v
}

echo ""
echo "========================================="
echo "Test suite completed!"
echo "========================================="
