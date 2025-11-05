"""Unit tests for fastq_random_sample.py"""

import unittest
import tempfile
import os
import gzip
import sys

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from fastq_random_sample import subsample_fastq, open_maybe_gzip


class TestFastqRandomSample(unittest.TestCase):
    """Test suite for fastq_random_sample.py"""

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

    def count_fastq_reads(self, fastq_path):
        """Helper: Count reads in a FASTQ file"""
        count = 0
        with open_maybe_gzip(fastq_path, "rt") as f:
            while True:
                header = f.readline()
                if not header:
                    break
                f.readline()  # sequence
                f.readline()  # plus
                f.readline()  # quality
                count += 1
        return count

    def test_subsample_basic(self):
        """Test basic subsampling functionality"""
        output_file = os.path.join(self.temp_dir, "sampled.fastq.gz")

        # Sample 10 reads
        subsample_fastq(self.test_fastq, output_file, num_reads=10, seed=42)

        # Verify output exists
        self.assertTrue(os.path.exists(output_file))

        # Verify correct number of reads
        self.assertEqual(self.count_fastq_reads(output_file), 10)

    def test_subsample_reproducible(self):
        """Test that same seed produces same results"""
        output1 = os.path.join(self.temp_dir, "sample1.fastq.gz")
        output2 = os.path.join(self.temp_dir, "sample2.fastq.gz")

        # Sample with same seed twice
        subsample_fastq(self.test_fastq, output1, num_reads=5, seed=123)
        subsample_fastq(self.test_fastq, output2, num_reads=5, seed=123)

        # Read both outputs
        with gzip.open(output1, 'rt') as f1, gzip.open(output2, 'rt') as f2:
            content1 = f1.read()
            content2 = f2.read()

        # Should be identical
        self.assertEqual(content1, content2)

    def test_subsample_different_seeds(self):
        """Test that different seeds produce different results"""
        output1 = os.path.join(self.temp_dir, "sample1.fastq.gz")
        output2 = os.path.join(self.temp_dir, "sample2.fastq.gz")

        # Sample with different seeds
        subsample_fastq(self.test_fastq, output1, num_reads=10, seed=1)
        subsample_fastq(self.test_fastq, output2, num_reads=10, seed=2)

        # Read both outputs
        with gzip.open(output1, 'rt') as f1, gzip.open(output2, 'rt') as f2:
            content1 = f1.read()
            content2 = f2.read()

        # Should be different (very high probability)
        self.assertNotEqual(content1, content2)

    def test_subsample_plain_fastq(self):
        """Test subsampling plain (non-gzipped) FASTQ"""
        # Create plain FASTQ from gzipped
        plain_input = os.path.join(self.temp_dir, "input.fastq")
        with gzip.open(self.test_fastq, 'rt') as fin, open(plain_input, 'w') as fout:
            # Copy first 20 reads
            for _ in range(20):
                for _ in range(4):
                    line = fin.readline()
                    if not line:
                        break
                    fout.write(line)

        output_file = os.path.join(self.temp_dir, "sampled.fastq")
        subsample_fastq(plain_input, output_file, num_reads=5, seed=42)

        self.assertTrue(os.path.exists(output_file))
        self.assertEqual(self.count_fastq_reads(output_file), 5)

    def test_subsample_too_many_reads_raises_error(self):
        """Test that requesting more reads than available raises ValueError"""
        output_file = os.path.join(self.temp_dir, "sampled.fastq.gz")

        # Count total reads in test file
        total_reads = self.count_fastq_reads(self.test_fastq)

        # Request more than available
        with self.assertRaises(ValueError) as context:
            subsample_fastq(self.test_fastq, output_file, num_reads=total_reads + 100, seed=42)

        self.assertIn("Requested", str(context.exception))
        self.assertIn("contains only", str(context.exception))

    def test_open_maybe_gzip_gzipped(self):
        """Test open_maybe_gzip correctly handles gzipped files"""
        with open_maybe_gzip(self.test_fastq, "rt") as f:
            first_line = f.readline()

        # Should start with @ (FASTQ header)
        self.assertTrue(first_line.startswith('@'))

    def test_open_maybe_gzip_plain(self):
        """Test open_maybe_gzip correctly handles plain files"""
        plain_file = os.path.join(self.temp_dir, "test.fastq")
        with open(plain_file, 'w') as f:
            f.write("@read1\nATCG\n+\nIIII\n")

        with open_maybe_gzip(plain_file, "rt") as f:
            first_line = f.readline()

        self.assertEqual(first_line.strip(), "@read1")


if __name__ == '__main__':
    unittest.main()
