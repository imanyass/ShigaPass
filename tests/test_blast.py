"""Tests for blast module"""

import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock
import shutil

from shigapass.blast import Blast


class TestBlast:
    @pytest.fixture
    def blast_instance(self):
        """Create a Blast instance for testing"""
        return Blast(db="/fake/db", threads=2)

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing"""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_blast_init(self):
        """Test Blast class initialization"""
        blast = Blast(db="/test/db", threads=4)
        assert blast.db == "/test/db"
        assert blast.threads == 4

    def test_get_hits(self, blast_instance, temp_dir):
        """Test get_hits method"""
        # Create test data
        sample_dir = os.path.join(temp_dir, "sample1")
        os.makedirs(sample_dir)

        test_records = "query1\tgene1_100\t99.0\t150\n"
        test_records += "query2\tgene2_50\t98.5\t120\n"
        test_records += "query3\tgene1_75\t97.0\t130\n"

        records_file = os.path.join(sample_dir, "test_allrecords.txt")
        with open(records_file, "w") as f:
            f.write(test_records)

        # Test the method
        result = blast_instance.get_hits(temp_dir, "sample1", "test")

        # Check results
        assert "gene1" in result
        assert "gene2" in result
        assert result["gene1"] == 2  # gene1 appears twice
        assert result["gene2"] == 1  # gene2 appears once

        # Check output file was created
        hits_file = os.path.join(sample_dir, "test_hits.txt")
        assert os.path.exists(hits_file)

    def test_filter_blast(self, blast_instance, temp_dir):
        """Test filter_blast method"""
        # Create test BLAST output
        blast_output = "query1\tsubject1_150\t99.0\t150\t1\t150\t1\t150\t1e-50\t200\n"
        blast_output += "query2\tsubject2_100\t95.0\t80\t1\t80\t1\t80\t1e-30\t150\n"
        blast_output += "query3\tsubject3_200\t97.0\t90\t1\t90\t1\t90\t1e-40\t180\n"

        blast_file = os.path.join(temp_dir, "blast_output.txt")
        filtered_file = os.path.join(temp_dir, "filtered_output.txt")

        with open(blast_file, "w") as f:
            f.write(blast_output)

        # Test filtering
        blast_instance.filter_blast(blast_file, filtered_file, identity=98, coverage=95)

        # Check filtered output
        with open(filtered_file, "r") as f:
            lines = f.readlines()

        # Only first line should pass (99% identity, 100% coverage)
        assert len(lines) == 1
        assert "subject1_150" in lines[0]

    @patch('subprocess.run')
    def test_run_blastn(self, mock_subprocess, blast_instance, temp_dir):
        """Test run_blastn method"""
        # Create sample directory
        sample_dir = os.path.join(temp_dir, "sample1")
        os.makedirs(sample_dir)

        # Mock successful subprocess run
        mock_subprocess.return_value = MagicMock()

        # Create empty blast output file (would normally be created by blastn)
        blast_output = os.path.join(sample_dir, "test_blastout.txt")
        with open(blast_output, "w") as f:
            f.write("query1\tsubject1_150\t99.0\t150\t1\t150\t1\t150\t1e-50\t200\n")

        with patch.object(blast_instance, 'filter_blast') as mock_filter:
            blast_instance.run_blastn(
                temp_dir, "sample1", "/fake/input.fasta", 
                "db.fasta", "test", 98, 95
            )

            # Check that subprocess was called with correct arguments
            mock_subprocess.assert_called_once()
            args = mock_subprocess.call_args[0][0]
            assert "blastn" in args
            assert "/fake/input.fasta" in args

            # Check that filter_blast was called
            mock_filter.assert_called_once()

    @patch('subprocess.run')
    def test_make_blast_db(self, mock_subprocess, blast_instance, temp_dir):
        """Test make_blast_db method"""
        # Create fake database structure
        db_dir = os.path.join(temp_dir, "IPAH")
        os.makedirs(db_dir)

        fasta_file = os.path.join(db_dir, "test.fasta")
        with open(fasta_file, "w") as f:
            f.write(">seq1\nATCG\n")

        # Set blast instance db to temp directory
        blast_instance.db = temp_dir

        # Mock successful subprocess run
        mock_subprocess.return_value = MagicMock()

        # Test make_blast_db
        blast_instance.make_blast_db()

        # Check that makeblastdb was called
        mock_subprocess.assert_called()
        args = mock_subprocess.call_args[0][0]
        assert "makeblastdb" in args
        assert "-dbtype" in args
        assert "nucl" in args
