"""Unit tests for map_fastq_to_fasta.py"""

import unittest
import tempfile
import os
import sys

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from map_fastq_to_fasta import (
    is_fastq_or_gz,
    fastq_stem,
    is_fasta_like,
    extract_species,
    parse_cigar,
    get_optional_tag,
    is_unmapped,
    is_secondary_or_supp
)


class TestMapFastqToFasta(unittest.TestCase):
    """Test suite for map_fastq_to_fasta.py"""

    def test_is_fastq_or_gz(self):
        """Test FASTQ file extension detection"""
        self.assertTrue(is_fastq_or_gz("reads.fastq"))
        self.assertTrue(is_fastq_or_gz("reads.fastq.gz"))
        self.assertTrue(is_fastq_or_gz("reads.fq"))
        self.assertTrue(is_fastq_or_gz("reads.fq.gz"))
        self.assertFalse(is_fastq_or_gz("reads.fasta"))
        self.assertFalse(is_fastq_or_gz("reads.txt"))
        self.assertFalse(is_fastq_or_gz("reads"))

    def test_fastq_stem(self):
        """Test FASTQ filename stem extraction"""
        self.assertEqual(fastq_stem("sample.fastq"), "sample")
        self.assertEqual(fastq_stem("sample.fastq.gz"), "sample")
        self.assertEqual(fastq_stem("sample.fq"), "sample")
        self.assertEqual(fastq_stem("sample.fq.gz"), "sample")
        self.assertEqual(fastq_stem("my.sample.fastq.gz"), "my.sample")
        self.assertEqual(fastq_stem("sample.txt"), "sample.txt")  # Not a FASTQ

    def test_is_fasta_like(self):
        """Test FASTA file extension detection"""
        self.assertTrue(is_fasta_like("ref.fasta"))
        self.assertTrue(is_fasta_like("ref.fa"))
        self.assertTrue(is_fasta_like("ref.fna"))
        self.assertFalse(is_fasta_like("ref.fastq"))
        self.assertFalse(is_fasta_like("ref.txt"))
        self.assertFalse(is_fasta_like("ref"))

    def test_extract_species(self):
        """Test species name extraction from FASTA headers"""
        # Valid species name
        self.assertEqual(
            extract_species(">NC_001234 Homo sapiens chromosome 1"),
            "Homo sapiens"
        )
        self.assertEqual(
            extract_species(">ref|NC_001234 Escherichia coli strain K12"),
            "Escherichia coli"
        )

        # No species name
        self.assertIsNone(extract_species(">sequence_id"))
        self.assertIsNone(extract_species(">123456"))

    def test_parse_cigar_basic(self):
        """Test CIGAR string parsing - basic operations"""
        # Match
        result = parse_cigar("100M")
        self.assertEqual(result['M'], 100)
        self.assertEqual(result['I'], 0)
        self.assertEqual(result['D'], 0)

        # Match + Insertion
        result = parse_cigar("50M10I40M")
        self.assertEqual(result['M'], 90)
        self.assertEqual(result['I'], 10)
        self.assertEqual(result['D'], 0)

        # Match + Deletion
        result = parse_cigar("30M5D70M")
        self.assertEqual(result['M'], 100)
        self.assertEqual(result['D'], 5)

        # Complex
        result = parse_cigar("10M2I5M3D20M")
        self.assertEqual(result['M'], 35)
        self.assertEqual(result['I'], 2)
        self.assertEqual(result['D'], 3)

    def test_parse_cigar_special_cases(self):
        """Test CIGAR string parsing - special cases"""
        # Empty/missing CIGAR
        result = parse_cigar("*")
        self.assertEqual(result['M'], 0)

        result = parse_cigar("")
        self.assertEqual(result['M'], 0)

        # Soft clipping
        result = parse_cigar("5S50M5S")
        self.assertEqual(result['S'], 10)
        self.assertEqual(result['M'], 50)

        # Hard clipping
        result = parse_cigar("10H50M10H")
        self.assertEqual(result['H'], 20)
        self.assertEqual(result['M'], 50)

    def test_parse_cigar_all_operations(self):
        """Test CIGAR string parsing - all operation types"""
        cigar = "5M2I3M1D10M2N5M3S4H"
        result = parse_cigar(cigar)

        self.assertEqual(result['M'], 23)  # 5 + 3 + 10 + 5
        self.assertEqual(result['I'], 2)
        self.assertEqual(result['D'], 1)
        self.assertEqual(result['N'], 2)
        self.assertEqual(result['S'], 3)
        self.assertEqual(result['H'], 4)

    def test_get_optional_tag_integer(self):
        """Test SAM optional tag extraction - integer values"""
        tags = ["NM:i:5", "AS:i:100", "XS:i:50"]

        self.assertEqual(get_optional_tag(tags, "NM"), 5)
        self.assertEqual(get_optional_tag(tags, "AS"), 100)
        self.assertEqual(get_optional_tag(tags, "XS"), 50)
        self.assertIsNone(get_optional_tag(tags, "XX"))  # Not present

    def test_get_optional_tag_string(self):
        """Test SAM optional tag extraction - string values"""
        tags = ["RG:Z:sample1", "MD:Z:50A10T30"]

        self.assertEqual(get_optional_tag(tags, "RG"), "sample1")
        self.assertEqual(get_optional_tag(tags, "MD"), "50A10T30")

    def test_get_optional_tag_invalid(self):
        """Test SAM optional tag extraction - invalid formats"""
        tags = ["INVALID", "NM:5"]  # Missing type field

        self.assertIsNone(get_optional_tag(tags, "NM"))

    def test_get_optional_tag_invalid_integer(self):
        """Test SAM optional tag extraction - invalid integer"""
        tags = ["NM:i:abc"]  # Not a valid integer

        self.assertIsNone(get_optional_tag(tags, "NM"))

    def test_is_unmapped(self):
        """Test SAM flag unmapped detection"""
        # Unmapped (flag 4)
        self.assertTrue(is_unmapped(4))
        self.assertTrue(is_unmapped(4 | 1))  # Unmapped + paired
        self.assertTrue(is_unmapped(4 | 16))  # Unmapped + reverse

        # Mapped
        self.assertFalse(is_unmapped(0))
        self.assertFalse(is_unmapped(16))  # Reverse strand
        self.assertFalse(is_unmapped(99))  # Paired, first in pair

    def test_is_secondary_or_supp(self):
        """Test SAM flag secondary/supplementary detection"""
        # Secondary alignment (flag 256)
        self.assertTrue(is_secondary_or_supp(256))
        self.assertTrue(is_secondary_or_supp(256 | 16))  # Secondary + reverse

        # Supplementary alignment (flag 2048)
        self.assertTrue(is_secondary_or_supp(2048))
        self.assertTrue(is_secondary_or_supp(2048 | 16))  # Supplementary + reverse

        # Both
        self.assertTrue(is_secondary_or_supp(256 | 2048))

        # Primary alignment
        self.assertFalse(is_secondary_or_supp(0))
        self.assertFalse(is_secondary_or_supp(16))
        self.assertFalse(is_secondary_or_supp(99))


class TestMapFastqToFastaIntegration(unittest.TestCase):
    """Integration tests for map_fastq_to_fasta.py (requires minimap2)"""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures"""
        cls.test_data_dir = os.path.join(os.path.dirname(__file__), '..', 'test_data')

    def setUp(self):
        """Set up temporary directory for each test"""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files after each test"""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @unittest.skipUnless(
        os.system("which minimap2 > /dev/null 2>&1") == 0,
        "minimap2 not available"
    )
    def test_map_reads_minimap2(self):
        """Test full mapping pipeline with minimap2 (if available)"""
        from map_fastq_to_fasta import map_reads

        target_folder = os.path.join(self.test_data_dir, 'resources', 'panel')
        query_folder = os.path.join(self.test_data_dir, 'data')

        summary_file = map_reads(
            target_folder=target_folder,
            query_folder=query_folder,
            threads=2,
            software="minimap2",
            output_dir=self.temp_dir,
            min_identity=80.0,
            max_size_diff=20.0
        )

        # Verify summary file created
        self.assertTrue(os.path.exists(summary_file))

        # Verify CSV has content
        with open(summary_file) as f:
            lines = f.readlines()
        self.assertGreater(len(lines), 1)  # Header + at least one data row

        # Verify header
        self.assertIn("Query_File", lines[0])
        self.assertIn("Target_File", lines[0])


if __name__ == '__main__':
    unittest.main()
