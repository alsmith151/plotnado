"""
Tests for the GenomicRegion class.
"""

from plotnado.tracks import GenomicRegion


class TestGenomicRegion:
    """Test cases for GenomicRegion class."""

    def test_creation(self):
        """Test creating a GenomicRegion instance."""
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")
        assert gr.chromosome == "chr1"
        assert gr.start == 1000
        assert gr.end == 2000
        assert gr.strand == "+"

    def test_length_property(self):
        """Test the length property."""
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")
        assert gr.length == 1000

    def test_center_property(self):
        """Test the center property."""
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")
        assert gr.center == 1500

    def test_string_representation(self):
        """Test the string representation."""
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")
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
        gr_dict = {"chromosome": "chr1", "start": 1000, "end": 2000, "strand": "+"}
        gr = GenomicRegion.from_dict(gr_dict)
        assert gr.chromosome == "chr1"
        assert gr.start == 1000
        assert gr.end == 2000
        assert gr.strand == "+"

    def test_extend(self):
        """Test extending the region."""
        # Positive strand
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")
        extended = gr.extend(upstream=100, downstream=200)
        assert extended.start == 900
        assert extended.end == 2200
        assert extended.strand == "+"

        # Negative strand
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="-")
        extended = gr.extend(upstream=100, downstream=200)
        assert (
            extended.start == 800
        )  # downstream on negative strand means lower coordinates
        assert (
            extended.end == 2100
        )  # upstream on negative strand means higher coordinates
        assert extended.strand == "-"

        # Negative start prevention
        gr = GenomicRegion(chromosome="chr1", start=100, end=200, strand="+")
        extended = gr.extend(upstream=1000)
        assert extended.start == 0
