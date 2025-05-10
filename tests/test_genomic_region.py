"""
Tests for the GenomicRegion class.
"""

import pytest
from plotnado.tracks import GenomicRegion


class TestGenomicRegion:
    """Test cases for GenomicRegion class."""

    def test_creation(self):
        """Test creating a GenomicRegion instance."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        assert gr.chromosome == "chr1"
        assert gr.start == 1000
        assert gr.end == 2000
        assert gr.strand == "+"

    def test_length_property(self):
        """Test the length property."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        assert gr.length == 1000

    def test_center_property(self):
        """Test the center property."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        assert gr.center == 1500

    def test_string_representation(self):
        """Test the string representation."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        assert str(gr) == "chr1:1000-2000(+)"

    def test_from_str(self):
        """Test creating from a string representation."""
        gr_str = "chr1:1000-2000(+)"
        gr = GenomicRegion.from_str(gr_str)
        assert gr.chromosome == "chr1"
        assert gr.start == 1000
        assert gr.end == 2000
        assert gr.strand == "+"

    def test_from_tuple(self):
        """Test creating from a tuple."""
        gr_tuple = ("chr1", 1000, 2000, "+")
        gr = GenomicRegion.from_tuple(gr_tuple)
        assert gr.chromosome == "chr1"
        assert gr.start == 1000
        assert gr.end == 2000
        assert gr.strand == "+"

    def test_from_list(self):
        """Test creating from a list."""
        gr_list = ["chr1", 1000, 2000, "+"]
        gr = GenomicRegion.from_list(gr_list)
        assert gr.chromosome == "chr1"
        assert gr.start == 1000
        assert gr.end == 2000
        assert gr.strand == "+"

    def test_from_dict(self):
        """Test creating from a dictionary."""
        gr_dict = {
            "chromosome": "chr1",
            "start": 1000,
            "end": 2000,
            "strand": "+"
        }
        gr = GenomicRegion.from_dict(gr_dict)
        assert gr.chromosome == "chr1"
        assert gr.start == 1000
        assert gr.end == 2000
        assert gr.strand == "+"
