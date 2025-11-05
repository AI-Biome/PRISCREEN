"""Unit tests for fastq_error_simulation.py"""

import unittest
import tempfile
import os
import gzip
import sys

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from fastq_error_simulation import simulate_errors, mutate_base, open_maybe_gzip


class TestFastqErrorSimulation(unittest.TestCase):
    """Test suite for fastq_error_simulation.py"""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests"""
        cls.test_data_dir = os.path.join(os.path.dirname(__file__), '..', 'test_data', 'data')
        cls.test_fastq = os.path.join(cls.test_data_dir, 'Test1.fastq.gz')

        # Verify test data exists
        if not os.path.exists(cls.test_fastq):
            raise FileNotFoundError(f"Test data not found: {cls.test_fastq}")

    def setUp(self):
        """Set up temporary directory for each test"""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files after each test"""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def count_differences(self, seq1, seq2):
        """Helper: Count number of differences between two sequences"""
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

    def test_mutate_base_changes_base(self):
        """Test that mutate_base returns a different base"""
        original = "A"
        mutated = mutate_base(original)

        # Should be different
        self.assertNotEqual(original, mutated)
        # Should be a valid DNA base
        self.assertIn(mutated, ["C", "G", "T", "N"])

    def test_mutate_base_preserves_case(self):
        """Test that mutate_base preserves case"""
        # Uppercase
        mutated_upper = mutate_base("A")
        self.assertTrue(mutated_upper.isupper())

        # Lowercase
        mutated_lower = mutate_base("a")
        self.assertTrue(mutated_lower.islower())

    def test_simulate_errors_zero_rate(self):
        """Test that 0% error rate produces identical output"""
        output_file = os.path.join(self.temp_dir, "mutated.fastq.gz")

        simulate_errors(self.test_fastq, output_file, error_rate=0.0, seed=42)

        # Read original and output
        with gzip.open(self.test_fastq, 'rt') as f_orig, gzip.open(output_file, 'rt') as f_out:
            orig_content = f_orig.read()
            out_content = f_out.read()

        # Should be identical
        self.assertEqual(orig_content, out_content)

    def test_simulate_errors_reproducible(self):
        """Test that same seed produces same results"""
        output1 = os.path.join(self.temp_dir, "mutated1.fastq.gz")
        output2 = os.path.join(self.temp_dir, "mutated2.fastq.gz")

        # Simulate with same seed twice
        simulate_errors(self.test_fastq, output1, error_rate=0.01, seed=123)
        simulate_errors(self.test_fastq, output2, error_rate=0.01, seed=123)

        # Read both outputs
        with gzip.open(output1, 'rt') as f1, gzip.open(output2, 'rt') as f2:
            content1 = f1.read()
            content2 = f2.read()

        # Should be identical
        self.assertEqual(content1, content2)

    def test_simulate_errors_introduces_mutations(self):
        """Test that error rate > 0 introduces mutations"""
        output_file = os.path.join(self.temp_dir, "mutated.fastq.gz")

        simulate_errors(self.test_fastq, output_file, error_rate=0.05, seed=42)

        # Read first read from original and mutated
        with gzip.open(self.test_fastq, 'rt') as f_orig, gzip.open(output_file, 'rt') as f_out:
            # Skip to sequence line
            f_orig.readline()
            orig_seq = f_orig.readline().strip()
            f_out.readline()
            mut_seq = f_out.readline().strip()

        # Should have some differences (with 5% error rate over ~80bp, very likely)
        differences = self.count_differences(orig_seq, mut_seq)
        self.assertGreater(differences, 0)

    def test_simulate_errors_quality_updated(self):
        """Test that quality scores are updated for mutated bases"""
        output_file = os.path.join(self.temp_dir, "mutated.fastq.gz")

        simulate_errors(
            self.test_fastq,
            output_file,
            error_rate=0.1,  # High error rate
            seed=42,
            low_quality_char='!',
            preserve_quality=False
        )

        # Read quality line
        with gzip.open(output_file, 'rt') as f:
            f.readline()  # header
            seq = f.readline().strip()
            f.readline()  # plus
            qual = f.readline().strip()

        # Should have low quality characters where mutations occurred
        self.assertIn('!', qual)

    def test_simulate_errors_preserve_quality(self):
        """Test that preserve_quality keeps original quality scores"""
        output_file = os.path.join(self.temp_dir, "mutated.fastq.gz")

        simulate_errors(
            self.test_fastq,
            output_file,
            error_rate=0.1,
            seed=42,
            preserve_quality=True
        )

        # Read original and mutated quality
        with gzip.open(self.test_fastq, 'rt') as f_orig, gzip.open(output_file, 'rt') as f_out:
            for _ in range(3):
                f_orig.readline()
                f_out.readline()
            orig_qual = f_orig.readline().strip()
            out_qual = f_out.readline().strip()

        # Quality should be identical even if sequence changed
        self.assertEqual(orig_qual, out_qual)

    def test_simulate_errors_plain_fastq(self):
        """Test error simulation on plain (non-gzipped) FASTQ"""
        # Create plain FASTQ
        plain_input = os.path.join(self.temp_dir, "input.fastq")
        plain_output = os.path.join(self.temp_dir, "output.fastq")

        with gzip.open(self.test_fastq, 'rt') as fin, open(plain_input, 'w') as fout:
            # Copy first 10 reads
            for _ in range(10):
                for _ in range(4):
                    line = fin.readline()
                    if not line:
                        break
                    fout.write(line)

        simulate_errors(plain_input, plain_output, error_rate=0.01, seed=42)

        # Should produce output
        self.assertTrue(os.path.exists(plain_output))
        self.assertGreater(os.path.getsize(plain_output), 0)

    def test_simulate_errors_invalid_rate_raises_error(self):
        """Test that invalid error rates raise ValueError"""
        output_file = os.path.join(self.temp_dir, "mutated.fastq.gz")

        # Test negative rate
        with self.assertRaises(ValueError) as context:
            simulate_errors(self.test_fastq, output_file, error_rate=-0.1)
        self.assertIn("between 0 and 1", str(context.exception))

        # Test rate > 1
        with self.assertRaises(ValueError) as context:
            simulate_errors(self.test_fastq, output_file, error_rate=1.5)
        self.assertIn("between 0 and 1", str(context.exception))

    def test_simulate_errors_malformed_fastq_raises_error(self):
        """Test that malformed FASTQ raises ValueError"""
        malformed_file = os.path.join(self.temp_dir, "malformed.fastq")
        output_file = os.path.join(self.temp_dir, "output.fastq")

        # Create malformed FASTQ (missing + line)
        with open(malformed_file, 'w') as f:
            f.write("@read1\n")
            f.write("ATCG\n")
            f.write("IIII\n")  # Missing + line

        with self.assertRaises(ValueError) as context:
            simulate_errors(malformed_file, output_file, error_rate=0.01)
        self.assertIn("Malformed FASTQ", str(context.exception))


if __name__ == '__main__':
    unittest.main()
