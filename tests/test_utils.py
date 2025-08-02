"""Tests for utils module"""

import os
import tempfile
import pytest
import shutil
import pandas as pd

from shigapass.utils import dict_to_str, write_out, read_file, get_val_from_file, check_ascending


class TestUtils:
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing"""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_dict_to_str(self):
        """Test dict_to_str function"""
        test_dict = {"gene1": 10, "gene2": 5, "gene3": 15}
        result = dict_to_str(test_dict)

        lines = result.strip().split("\n")
        assert len(lines) == 3
        assert "gene1;10" in lines
        assert "gene2;5" in lines
        assert "gene3;15" in lines

    def test_write_out(self, temp_dir):
        """Test write_out function"""
        test_file = os.path.join(temp_dir, "test_output.txt")
        test_content = "This is test content\nLine 2\n"

        # Test writing
        write_out(test_content, test_file)

        # Verify file was created and content is correct
        assert os.path.exists(test_file)
        with open(test_file, "r") as f:
            content = f.read()
        assert content == test_content

        # Test appending
        additional_content = "Line 3\n"
        write_out(additional_content, test_file, "a")

        with open(test_file, "r") as f:
            content = f.read()
        assert content == test_content + additional_content

    def test_read_file(self, temp_dir):
        """Test read_file function"""
        test_file = os.path.join(temp_dir, "test_input.txt")
        test_content = "This is test content\nLine 2\nLine 3\n"

        # Create test file
        with open(test_file, "w") as f:
            f.write(test_content)

        # Test reading
        result = read_file(test_file)
        assert result == test_content

    def test_get_val_from_file(self, temp_dir):
        """Test get_val_from_file function"""
        test_file = os.path.join(temp_dir, "test.csv")
        test_content = "gene1;10;50\ngene2;20;30\ngene3;15;40\n"

        with open(test_file, "w") as f:
            f.write(test_content)

        # Test getting top value sorted by column 1 (hits)
        result = get_val_from_file(test_file, 1, 0, delimiter=";", ascending=False, num_of_rows=1)
        assert result == ["gene2"]  # gene2 has highest hits (20)

        # Test getting top 2 values
        result = get_val_from_file(test_file, 1, 0, delimiter=";", ascending=False, num_of_rows=2)
        assert result == ["gene2", "gene3"]  # gene2 (20), gene3 (15)

        # Test getting value from different column
        result = get_val_from_file(test_file, 1, 1, delimiter=";", ascending=False, num_of_rows=1)
        assert result == [20]  # The hits value for gene2

        # Test ascending sort
        result = get_val_from_file(test_file, 1, 0, delimiter=";", ascending=True, num_of_rows=1)
        assert result == ["gene1"]  # gene1 has lowest hits (10)

    def test_get_val_from_file_empty(self, temp_dir):
        """Test get_val_from_file with empty file"""
        test_file = os.path.join(temp_dir, "empty.csv")

        with open(test_file, "w") as f:
            f.write("")

        result = get_val_from_file(test_file, 1, 0)
        assert result == ""

    def test_check_ascending(self, temp_dir):
        """Test check_ascending function"""
        test_file = os.path.join(temp_dir, "test_blast.txt")

        # Test ascending order (column 8 < column 9)
        ascending_content = "query1\tsubject1\t99.0\t150\t1\t150\t1\t150\t100\t200\n"
        with open(test_file, "w") as f:
            f.write(ascending_content)

        result = check_ascending(test_file, 8, 9)
        assert result is True

        # Test descending order (column 8 > column 9)
        descending_content = "query1\tsubject1\t99.0\t150\t1\t150\t1\t150\t200\t100\n"
        with open(test_file, "w") as f:
            f.write(descending_content)

        result = check_ascending(test_file, 8, 9)
        assert result is False

        # Test equal values
        equal_content = "query1\tsubject1\t99.0\t150\t1\t150\t1\t150\t150\t150\n"
        with open(test_file, "w") as f:
            f.write(equal_content)

        result = check_ascending(test_file, 8, 9)
        assert result is True
