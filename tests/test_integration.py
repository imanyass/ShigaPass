"""Integration tests for ShigaPass"""

import os
import tempfile
import pytest
import shutil
from unittest.mock import patch, MagicMock

from shigapass.main import OptionsParser
from shigapass.blast import Blast
from shigapass.typing import Typing


class TestIntegration:
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing"""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_options_parser_init(self):
        """Test OptionsParser initialization"""
        parser = OptionsParser("1.5.0")
        assert parser.version == "1.5.0"

    @patch('subprocess.run')
    def test_full_pipeline_negative_ipah(self, mock_subprocess, temp_dir):
        """Test full pipeline with negative ipaH result"""
        # Create test input files
        input_file = os.path.join(temp_dir, "test.fasta")
        with open(input_file, "w") as f:
            f.write(">test_sequence\nATCGATCGATCG\n")

        list_file = os.path.join(temp_dir, "input_list.txt")
        with open(list_file, "w") as f:
            f.write(input_file + "\n")

        output_dir = os.path.join(temp_dir, "output")
        db_dir = os.path.join(temp_dir, "db")

        # Create mock database structure
        os.makedirs(os.path.join(db_dir, "IPAH"))
        os.makedirs(os.path.join(db_dir, "RFB"))

        # Create RFB_hits_count.csv
        rfb_hits_file = os.path.join(db_dir, "RFB_hits_count.csv")
        with open(rfb_hits_file, "w") as f:
            f.write("A1;50\nB1;40\nC1;30\n")

        # Mock subprocess to avoid actual BLAST calls
        mock_subprocess.return_value = MagicMock()

        # Create Blast and Typing instances
        blast = Blast(db_dir, 2)
        typing = Typing(blast)

        # Mock blast results to simulate negative ipaH
        def mock_search(outdir, input_filename, input_fasta_file, search_parameters):
            # Create empty allrecords file
            sample_dir = os.path.join(outdir, input_filename)
            marker = search_parameters[0]["db_marker"]
            allrecords_file = os.path.join(sample_dir, f"{marker}_allrecords.txt")
            with open(allrecords_file, "w") as f:
                f.write("")  # Empty file = no hits
            return {}, 0

        blast.search = mock_search

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Run the pipeline
        typing.run(list_file, output_dir, keep_files=True)

        # Check that summary file was created
        summary_file = os.path.join(output_dir, "ShigaPass_summary.csv")
        assert os.path.exists(summary_file)

        # Check summary content
        with open(summary_file, "r") as f:
            lines = f.readlines()

        assert len(lines) >= 2  # Header + at least one result
        assert "Not Shigella/EIEC" in lines[1]

    @patch('subprocess.run')
    def test_full_pipeline_positive_ipah(self, mock_subprocess, temp_dir):
        """Test full pipeline with positive ipaH result"""
        # Create test input files
        input_file = os.path.join(temp_dir, "test.fasta")
        with open(input_file, "w") as f:
            f.write(">test_sequence\nATCGATCGATCG\n")

        list_file = os.path.join(temp_dir, "input_list.txt")
        with open(list_file, "w") as f:
            f.write(input_file + "\n")

        output_dir = os.path.join(temp_dir, "output")
        db_dir = os.path.join(temp_dir, "db")

        # Create mock database structure
        os.makedirs(os.path.join(db_dir, "IPAH"))
        os.makedirs(os.path.join(db_dir, "RFB"))
        os.makedirs(os.path.join(db_dir, "MLST"))
        os.makedirs(os.path.join(db_dir, "FLIC"))
        os.makedirs(os.path.join(db_dir, "CRISPR"))

        # Create required database files
        rfb_hits_file = os.path.join(db_dir, "RFB_hits_count.csv")
        with open(rfb_hits_file, "w") as f:
            f.write("A1;50\nB1;40\nC1;30\n")

        mlst_profiles_file = os.path.join(db_dir, "MLST", "ST_profiles.txt")
        with open(mlst_profiles_file, "w") as f:
            f.write("ST\tadk\tfumC\tgyrB\ticd\tmdh\tpurA\trecA\n")
            f.write("1\t1\t1\t1\t1\t1\t1\t1\n")

        # Mock subprocess to avoid actual BLAST calls
        mock_subprocess.return_value = MagicMock()

        # Create Blast and Typing instances
        blast = Blast(db_dir, 2)
        typing = Typing(blast)

        # Mock blast results to simulate positive results
        def mock_search(outdir, input_filename, input_fasta_file, search_parameters):
            sample_dir = os.path.join(outdir, input_filename)
            marker = search_parameters[0]["db_marker"]
            allrecords_file = os.path.join(sample_dir, f"{marker}_allrecords.txt")
            hits_file = os.path.join(sample_dir, f"{marker}_hits.txt")
            coverage_file = os.path.join(sample_dir, f"{marker}_hitscoverage.txt")

            if marker == "ipaH":
                # Positive ipaH result
                with open(allrecords_file, "w") as f:
                    f.write("query1\tipaH_150\t99.0\t150\n")
                with open(hits_file, "w") as f:
                    f.write("ipaH;10\n")
                with open(coverage_file, "w") as f:
                    f.write("ipaH;10;50;20.0\n")
                return {"ipaH": 10}, 20.0
            elif marker == "rfb":
                # A1 rfb result
                with open(allrecords_file, "w") as f:
                    f.write("query1\tA1_150\t99.0\t150\n")
                with open(hits_file, "w") as f:
                    f.write("A1;15\n")
                with open(coverage_file, "w") as f:
                    f.write("A1;15;50;30.0\n")
                return {"A1": 15}, 30.0
            else:
                # Empty results for other markers
                with open(allrecords_file, "w") as f:
                    f.write("")
                return {}, 0

        # Mock run_blastn for MLST
        def mock_run_blastn(outdir, input_filename, input_fasta_file, db_fasta, db_marker, identity, coverage):
            sample_dir = os.path.join(outdir, input_filename)
            allrecords_file = os.path.join(sample_dir, f"{db_marker}_allrecords.txt")
            with open(allrecords_file, "w") as f:
                if "adk" in db_marker:
                    f.write("query1\tadk-1:full\t100.0\t400\n")
                elif "fumC" in db_marker:
                    f.write("query1\tfumC-1:full\t100.0\t400\n")
                else:
                    f.write("query1\t{}-1:full\t100.0\t400\n".format(db_marker))

        # Mock get_val_from_file for flic determination
        def mock_get_val_from_file(csv_file, sort_column, value_column, delimiter=";", ascending=True, num_of_rows=1):
            if "flic" in csv_file:
                return ["fliC_H7"]  # Return mock flic result
            return [""]

        blast.search = mock_search
        blast.run_blastn = mock_run_blastn

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Mock utility functions and typing methods to avoid file dependencies
        with patch('shigapass.utils.get_val_from_file', side_effect=mock_get_val_from_file):
            with patch('shigapass.utils.check_ascending', return_value=True):
                with patch.object(typing, 'determine_flic', return_value="H7"):
                    with patch.object(typing, 'determine_crispr', return_value="1"):
                        # Run the pipeline
                        typing.run(list_file, output_dir, keep_files=True)

        # Check that summary file was created
        summary_file = os.path.join(output_dir, "ShigaPass_summary.csv")
        assert os.path.exists(summary_file)

        # Check summary content
        with open(summary_file, "r") as f:
            lines = f.readlines()

        assert len(lines) >= 2  # Header + at least one result
        result_line = lines[1]
        assert "A1" in result_line  # Should detect A1 serotype
        assert "ipaH+" in result_line  # Should be positive for ipaH

    def test_command_line_argument_parsing(self):
        """Test that command line arguments are properly parsed"""
        # This would typically test the CLI module
        from shigapass.cli import get_args

        parser = get_args()

        # Test that all required arguments are present
        args = parser.parse_args(["-l", "test.txt", "-o", "output", "-d", "db"])
        assert args.list_file == "test.txt"
        assert args.outdir == "output"
        assert args.db == "db"
        assert args.threads == 2
        assert not args.mkdb
        assert not args.keep_files
