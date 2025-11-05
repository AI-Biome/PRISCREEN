#!/usr/bin/env python3
"""
Integration tests for the Snakemake workflow
"""

import os
import sys
import subprocess
import shutil
import unittest
from pathlib import Path


class TestSnakemakeWorkflow(unittest.TestCase):
    """Test the complete Snakemake workflow"""

    @classmethod
    def setUpClass(cls):
        """Set up test environment"""
        cls.project_root = Path(__file__).parent.parent
        cls.test_data_dir = cls.project_root / "test_data"
        cls.test_output_dir = cls.project_root / "test_results"
        cls.snakefile = cls.project_root / "Snakefile"
        cls.test_config = cls.test_data_dir / "config" / "config.yaml"
        cls.test_samples = cls.test_data_dir / "config" / "samples.tsv"

        # Check required files exist
        if not cls.test_data_dir.exists():
            raise FileNotFoundError(
                f"Test data directory not found: {cls.test_data_dir}"
            )
        if not cls.test_config.exists():
            raise FileNotFoundError(
                f"Test config not found: {cls.test_config}"
            )
        if not cls.test_samples.exists():
            raise FileNotFoundError(
                f"Test samples file not found: {cls.test_samples}"
            )

    def setUp(self):
        """Clean up before each test"""
        if self.test_output_dir.exists():
            shutil.rmtree(self.test_output_dir)
        self.test_output_dir.mkdir(exist_ok=True)

    def tearDown(self):
        """Clean up after each test"""
        # Optionally keep test results for inspection
        # Uncomment to auto-clean:
        # if self.test_output_dir.exists():
        #     shutil.rmtree(self.test_output_dir)
        pass

    def run_snakemake(self, targets=None, extra_args=None, expect_success=True):
        """
        Helper function to run Snakemake commands

        Args:
            targets: List of target rules/files, or None for default target
            extra_args: Additional command line arguments
            expect_success: Whether to expect the command to succeed

        Returns:
            subprocess.CompletedProcess object
        """
        cmd = [
            "snakemake",
            "--configfile", str(self.test_config),
            "--config", f"samples_file={self.test_samples}",
            "--directory", str(self.project_root),
            "--cores", "2",
            "--use-conda",
        ]

        if targets:
            cmd.extend(targets)

        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(
            cmd,
            cwd=self.project_root,
            capture_output=True,
            text=True
        )

        if expect_success and result.returncode != 0:
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            self.fail(f"Snakemake command failed with return code {result.returncode}")

        return result

    def test_snakefile_syntax(self):
        """Test that the Snakefile has valid syntax"""
        result = self.run_snakemake(
            extra_args=["--dryrun", "--quiet"],
            expect_success=True
        )
        self.assertEqual(result.returncode, 0, "Snakefile syntax check failed")

    def test_workflow_dag(self):
        """Test that the workflow DAG can be generated"""
        result = self.run_snakemake(
            extra_args=["--dag"],
            expect_success=True
        )
        self.assertEqual(result.returncode, 0, "DAG generation failed")
        self.assertIn("digraph", result.stdout, "DAG output missing")

    def test_list_rules(self):
        """Test that rules can be listed"""
        result = subprocess.run(
            ["snakemake", "--list", "--configfile", str(self.test_config)],
            cwd=self.project_root,
            capture_output=True,
            text=True
        )
        self.assertEqual(result.returncode, 0, "List rules failed")

        # Check that expected rules are present
        expected_rules = [
            "all",
            "mmseqs_createdb",
            "map_to_panel",
            "amplicon_readnames",
            "amplicon_fastq",
            "cluster_vsearch",
            "consensus_polish",
            "mmseqs_search",
            "summarize_hits",
            "sample_summary"
        ]

        for rule in expected_rules:
            self.assertIn(rule, result.stdout, f"Rule '{rule}' not found in workflow")

    def test_dryrun_complete_workflow(self):
        """Test dry-run of the complete workflow"""
        result = self.run_snakemake(
            extra_args=["--dryrun"],
            expect_success=True
        )
        self.assertEqual(result.returncode, 0, "Dry-run failed")

    def test_rule_mmseqs_createdb(self):
        """Test the mmseqs_createdb rule in isolation"""
        target = "results/mmseqs/panelDB"
        result = self.run_snakemake(
            targets=[target],
            extra_args=["--forcerun", "mmseqs_createdb"],
            expect_success=True
        )
        self.assertEqual(result.returncode, 0, "mmseqs_createdb rule failed")
        self.assertTrue(
            (self.project_root / target).exists(),
            "mmseqs database not created"
        )

    def test_rule_map_to_panel(self):
        """Test the map_to_panel rule"""
        target = "results/map/Test1.panel.sorted.bam"
        result = self.run_snakemake(
            targets=[target],
            expect_success=True
        )
        self.assertEqual(result.returncode, 0, "map_to_panel rule failed")
        # Note: temp files might be deleted, so we check the workflow succeeded

    def test_complete_workflow_single_sample(self):
        """Test the complete workflow for a single sample"""
        target = "results/identify/Test1/species_summary.tsv"
        result = self.run_snakemake(
            targets=[target],
            expect_success=True
        )
        self.assertEqual(result.returncode, 0, "Complete workflow failed")

        output_file = self.project_root / target
        self.assertTrue(
            output_file.exists(),
            f"Expected output file not found: {output_file}"
        )

        # Check that output file has content
        self.assertGreater(
            output_file.stat().st_size,
            0,
            "Output file is empty"
        )

    def test_workflow_output_structure(self):
        """Test that the workflow produces expected output structure"""
        # Run workflow for one sample
        self.run_snakemake(
            targets=["results/identify/Test1/species_summary.tsv"],
            expect_success=True
        )

        # Check expected directories exist
        expected_dirs = [
            "results/mmseqs",
            "results/identify/Test1",
        ]

        for dir_path in expected_dirs:
            full_path = self.project_root / dir_path
            self.assertTrue(
                full_path.exists() and full_path.is_dir(),
                f"Expected directory not found: {dir_path}"
            )

    def test_workflow_idempotence(self):
        """Test that re-running the workflow doesn't redo completed work"""
        target = "results/identify/Test1/species_summary.tsv"

        # First run
        self.run_snakemake(targets=[target], expect_success=True)

        # Second run - should show nothing to do
        result = self.run_snakemake(targets=[target], expect_success=True)

        # Check that no rules were executed (idempotent)
        self.assertIn(
            "Nothing to be done",
            result.stdout,
            "Workflow should be idempotent but rules were re-executed"
        )

    def test_included_rules_accessible(self):
        """Test that all included rule files are accessible"""
        rule_files = [
            "rules/mapping.smk",
            "rules/clustering.smk",
            "rules/consensus.smk",
            "rules/identification.smk",
        ]

        for rule_file in rule_files:
            full_path = self.project_root / rule_file
            self.assertTrue(
                full_path.exists(),
                f"Rule file not found: {rule_file}"
            )


class TestSnakemakeWorkflowConfig(unittest.TestCase):
    """Test workflow configuration and parameters"""

    @classmethod
    def setUpClass(cls):
        cls.project_root = Path(__file__).parent.parent
        cls.test_config = cls.project_root / "test_data" / "config" / "config.yaml"

    def test_config_file_valid(self):
        """Test that config file is valid YAML"""
        import yaml

        with open(self.test_config) as f:
            config = yaml.safe_load(f)

        self.assertIsInstance(config, dict, "Config is not a dictionary")

        # Check required sections
        required_sections = [
            "multi_species",
            "cluster",
            "consensus",
            "identify",
            "threads"
        ]

        for section in required_sections:
            self.assertIn(
                section,
                config,
                f"Required config section missing: {section}"
            )

    def test_samples_file_valid(self):
        """Test that samples file is valid"""
        import pandas as pd

        samples_file = self.project_root / "test_data" / "config" / "samples.tsv"
        df = pd.read_csv(samples_file, sep="\t")

        # Check required columns
        self.assertIn("sample", df.columns, "Missing 'sample' column")
        self.assertIn("fastq", df.columns, "Missing 'fastq' column")

        # Check that files exist (relative to test_data)
        for _, row in df.iterrows():
            fastq_path = self.project_root / "test_data" / row["fastq"]
            self.assertTrue(
                fastq_path.exists(),
                f"FASTQ file not found: {fastq_path}"
            )


def suite():
    """Create test suite"""
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(TestSnakemakeWorkflow))
    test_suite.addTest(unittest.makeSuite(TestSnakemakeWorkflowConfig))
    return test_suite


if __name__ == "__main__":
    # Check if snakemake is available
    try:
        subprocess.run(
            ["snakemake", "--version"],
            capture_output=True,
            check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: snakemake not found. Please install snakemake to run these tests.")
        sys.exit(1)

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite())
    sys.exit(0 if result.wasSuccessful() else 1)
