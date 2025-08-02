"""Tests for typing module"""

import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock
import shutil
import pandas as pd

from shigapass.typing import Typing
from shigapass.blast import Blast


class TestTyping:
    @pytest.fixture
    def blast_mock(self):
        """Create a mock Blast instance"""
        blast = MagicMock(spec=Blast)
        blast.db = "/fake/db"
        return blast

    @pytest.fixture
    def typing_instance(self, blast_mock):
        """Create a Typing instance for testing"""
        return Typing(blast_mock)

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing"""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_typing_init(self, blast_mock):
        """Test Typing class initialization"""
        typing = Typing(blast_mock)
        assert typing.blast == blast_mock

    def test_get_db_array(self, typing_instance):
        """Test get_db_array method"""
        result = typing_instance.get_db_array("ipaH")
        assert len(result) == 1
        assert result[0]["db_marker"] == "ipaH"

        result = typing_instance.get_db_array("rfb")
        assert len(result) == 1
        assert result[0]["db_marker"] == "rfb"

        # Test with rfb filter
        result = typing_instance.get_db_array("additionalrfb", "C1")
        assert len(result) == 1
        assert result[0]["rfb"] == "C1"

    def test_determine_ipah_positive(self, typing_instance, temp_dir):
        """Test determine_ipah method with positive result"""
        # Mock blast.search to return positive result
        typing_instance.blast.search.return_value = ({"ipaH": 10}, 95.5)

        result = typing_instance.determine_ipah(temp_dir, "sample1", "/fake/input.fasta")

        assert result[0] == "ipaH+"
        assert result[1] == 10
        assert result[2] == 95.5

    def test_determine_ipah_negative(self, typing_instance, temp_dir):
        """Test determine_ipah method with negative result"""
        # Mock blast.search to return empty result
        typing_instance.blast.search.return_value = ({}, 0)

        result = typing_instance.determine_ipah(temp_dir, "sample1", "/fake/input.fasta")

        assert result[0] == "ipaH-"
        assert result[1] == 0
        assert result[2] == 0

    def test_determine_rfb_basic(self, typing_instance, temp_dir):
        """Test determine_rfb method with basic result"""
        # Mock blast.search to return A1 result
        typing_instance.blast.search.return_value = ({"A1": 15}, 90.0)

        with patch.object(typing_instance, 'determine_alt_rfb') as mock_alt_rfb:
            mock_alt_rfb.return_value = ("A1", 15, 90.0)
            with patch.object(typing_instance, 'check_multiple_rfbs') as mock_check:
                mock_check.return_value = ""

                result = typing_instance.determine_rfb(temp_dir, "sample1", "/fake/input.fasta")

                assert result[0] == "A1"  # rfb
                assert result[1] == 15    # rfb_hits
                assert result[2] == 90.0  # rfb_coverage
                assert result[3] == ""    # flexserotype
                assert result[4] == ""    # comments

    def test_determine_rfb_flexneri(self, typing_instance, temp_dir):
        """Test determine_rfb method with S. flexneri result"""
        # Mock blast.search to return B1-5 result
        typing_instance.blast.search.return_value = ({"B1-5": 20}, 85.0)

        with patch.object(typing_instance, 'determine_alt_rfb') as mock_alt_rfb:
            mock_alt_rfb.return_value = ("B1-5", 20, 85.0)
            with patch.object(typing_instance, 'get_poac') as mock_poac:
                mock_poac.return_value = "1a"
                with patch.object(typing_instance, 'check_multiple_rfbs') as mock_check:
                    mock_check.return_value = ""

                    result = typing_instance.determine_rfb(temp_dir, "sample1", "/fake/input.fasta")

                    assert result[0] == "B1-5"  # rfb
                    assert result[3] == "1a"    # flexserotype

    def test_mlst_search(self, typing_instance, temp_dir):
        """Test mlst_search method"""
        # Create mock MLST database
        mlst_data = {
            'ST': [1, 2, 3],
            'adk': [1, 2, 3],
            'fumC': [1, 2, 3],
            'gyrB': [1, 2, 3],
            'icd': [1, 2, 3],
            'mdh': [1, 2, 3],
            'purA': [1, 2, 3],
            'recA': [1, 2, 3]
        }

        with patch('pandas.read_csv') as mock_read_csv:
            mock_read_csv.return_value = pd.DataFrame(mlst_data)

            # Test exact match
            search_dict = {'adk': '1', 'fumC': '1', 'gyrB': '1', 'icd': '1', 'mdh': '1', 'purA': '1', 'recA': '1'}
            result = typing_instance.mlst_search(search_dict)
            assert result == 1

            # Test no match
            search_dict = {'adk': '99', 'fumC': '99', 'gyrB': '99', 'icd': '99', 'mdh': '99', 'purA': '99', 'recA': '99'}
            result = typing_instance.mlst_search(search_dict)
            assert result == "~"

    def test_determine_flic(self, typing_instance, temp_dir):
        """Test determine_flic method"""
        # Create sample directory and files
        sample_dir = os.path.join(temp_dir, "sample1")
        os.makedirs(sample_dir)

        # Mock blast.search
        typing_instance.blast.search.return_value = ({}, 0)

        # Create mock flic records file
        flic_records = "query1\tfliC_H7\t99.0\t150\t1\t150\t1\t150\t1e-50\t200\t99.0\t300\n"
        flic_file = os.path.join(sample_dir, "flic_allrecords.txt")
        with open(flic_file, "w") as f:
            f.write(flic_records)

        with patch('shigapass.utils.get_val_from_file') as mock_get_val:
            mock_get_val.return_value = ["fliC_H7"]

            result = typing_instance.determine_flic(temp_dir, "sample1", "/fake/input.fasta")
            assert result == "H7"

    def test_determine_crispr(self, typing_instance, temp_dir):
        """Test determine_crispr method"""
        # Create sample directory
        sample_dir = os.path.join(temp_dir, "sample1")
        os.makedirs(sample_dir)

        # Mock blast.search
        typing_instance.blast.search.return_value = ({}, 0)

        # Create mock crispr records file
        crispr_records = "query1\tCRISPR_1\t100.0\t30\t1\t30\t1\t30\t1e-10\t60\t100.0\t30\n"
        crispr_file = os.path.join(sample_dir, "crispr_allrecords.txt")
        with open(crispr_file, "w") as f:
            f.write(crispr_records)

        with patch('shigapass.utils.check_ascending') as mock_check:
            mock_check.return_value = True
            with patch('shigapass.utils.get_val_from_file') as mock_get_val:
                mock_get_val.return_value = ["CRISPR_1"]

                result = typing_instance.determine_crispr(temp_dir, "sample1", "/fake/input.fasta")
                assert result == "1"

    def test_get_poac(self, typing_instance, temp_dir):
        """Test get_poac method"""
        # Create sample directory and POAC hits file
        sample_dir = os.path.join(temp_dir, "sample1")
        os.makedirs(sample_dir)

        poac_hits = "gtrI;10\ngtrII;5\noac;3\n"
        poac_file = os.path.join(sample_dir, "POAC_hits.txt")
        with open(poac_file, "w") as f:
            f.write(poac_hits)

        result = typing_instance.get_poac(temp_dir, "sample1")

        # Check that score file was created
        score_file = os.path.join(sample_dir, "score.txt")
        assert os.path.exists(score_file)

        # Check flex summary was created
        flex_summary_file = os.path.join(temp_dir, "ShigaPass_Flex_summary.csv")
        assert os.path.exists(flex_summary_file)

        # The score should be 1 + 4 + 64 = 69 (no exact match in flex_score_array)
        assert result == "Unknown"

    def test_check_multiple_rfbs(self, typing_instance, temp_dir):
        """Test check_multiple_rfbs method"""
        # Create sample directory and coverage file
        sample_dir = os.path.join(temp_dir, "sample1")
        os.makedirs(sample_dir)

        # Test with multiple RFBs (3 lines)
        coverage_content = "A1;10;50;80.0\nB1;5;40;60.0\nC1;3;30;50.0\n"
        coverage_file = os.path.join(sample_dir, "rfb_hitscoverage.txt")
        with open(coverage_file, "w") as f:
            f.write(coverage_content)

        result = typing_instance.check_multiple_rfbs(temp_dir, "sample1", "rfb", "A1")
        # The method should count lines, not characters
        assert "More than one rfb is detected: 3" in result

        # Test with single RFB
        coverage_content = "A1;10;50;80.0\n"
        with open(coverage_file, "w") as f:
            f.write(coverage_content)

        result = typing_instance.check_multiple_rfbs(temp_dir, "sample1", "rfb", "A1")
        assert result == ""
