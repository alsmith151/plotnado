"""
Tests for the aesthetics classes.
"""

import pytest
from pydantic import ValidationError

from plotnado.tracks import (
    Aesthetics,
    BigwigAesthetics,
    ScaleBarAesthetics,
    GenesAesthetics,
)


class TestAesthetics:
    """Test cases for the Aesthetics class and factory methods."""

    def test_bigwig_factory(self):
        """Test the bigwig factory method."""
        aesthetics = Aesthetics.bigwig(
            color="red", fill=True, alpha=0.5, style="scatter", scatter_point_size=3.0
        )

        assert isinstance(aesthetics, BigwigAesthetics)
        assert aesthetics.color == "red"
        assert aesthetics.fill is True
        assert aesthetics.alpha == 0.5
        assert aesthetics.style == "scatter"
        assert aesthetics.scatter_point_size == 3.0

    def test_scalebar_factory(self):
        """Test the scalebar factory method."""
        aesthetics = Aesthetics.scalebar(
            color="blue", position="center", scale_distance=500
        )

        assert isinstance(aesthetics, ScaleBarAesthetics)
        assert aesthetics.color == "blue"
        assert aesthetics.position == "center"
        assert aesthetics.scale_distance == 500

    def test_genes_factory(self):
        """Test the genes factory method."""
        aesthetics = Aesthetics.genes(
            color="green",
            display="expanded",
            minimum_gene_length=1000,
            max_number_of_rows=10,
        )

        assert isinstance(aesthetics, GenesAesthetics)
        assert aesthetics.color == "green"
        assert aesthetics.display == "expanded"
        assert aesthetics.minimum_gene_length == 1000
        assert aesthetics.max_number_of_rows == 10


class TestBigwigAesthetics:
    """Test cases for the BigwigAesthetics class."""

    def test_default_values(self):
        """Test default values."""
        aesthetics = BigwigAesthetics()

        assert aesthetics.style == "fill"  # Changed for genome browser style
        assert aesthetics.color == "#2171b5"  # Genome browser blue
        assert aesthetics.fill is True
        assert aesthetics.alpha == 0.85  # Slightly transparent
        assert aesthetics.scatter_point_size == 1.0
        assert aesthetics.linewidth == 1.0
        assert aesthetics.min_value is None
        assert aesthetics.max_value is None
        assert aesthetics.plot_title is True
        assert aesthetics.title_location == "left"
        assert aesthetics.title_height == 0.5
        assert aesthetics.plot_scale is True
        assert aesthetics.scale_location == "left"
        assert aesthetics.scale_height == 0.5

    def test_custom_values(self):
        """Test custom values."""
        aesthetics = BigwigAesthetics(
            style="scatter",
            color="red",
            fill=False,
            alpha=0.5,
            scatter_point_size=2.0,
            linewidth=2.0,
            min_value=0.0,
            max_value=1.0,
            plot_title=False,
            title_location="right",
            title_height=0.8,
            plot_scale=False,
            scale_location="right",
            scale_height=0.8,
        )

        assert aesthetics.style == "scatter"
        assert aesthetics.color == "red"
        assert aesthetics.fill is False
        assert aesthetics.alpha == 0.5
        assert aesthetics.scatter_point_size == 2.0
        assert aesthetics.linewidth == 2.0
        assert aesthetics.min_value == 0.0
        assert aesthetics.max_value == 1.0
        assert aesthetics.plot_title is False
        assert aesthetics.title_location == "right"
        assert aesthetics.title_height == 0.8
        assert aesthetics.plot_scale is False
        assert aesthetics.scale_location == "right"
        assert aesthetics.scale_height == 0.8

    def test_invalid_style(self):
        """Test invalid style value."""
        with pytest.raises(ValidationError):
            BigwigAesthetics(style="invalid_style")

    def test_invalid_title_location(self):
        """Test invalid title_location value."""
        with pytest.raises(ValidationError):
            BigwigAesthetics(title_location="invalid_location")

    def test_invalid_scale_location(self):
        """Test invalid scale_location value."""
        with pytest.raises(ValidationError):
            BigwigAesthetics(scale_location="invalid_location")


class TestScaleBarAesthetics:
    """Test cases for the ScaleBarAesthetics class."""

    def test_default_values(self):
        """Test default values."""
        aesthetics = ScaleBarAesthetics()

        assert aesthetics.style == "std"
        assert aesthetics.color == "#333333"  # Dark gray
        assert aesthetics.position == "left"
        assert aesthetics.scale_distance is None
        assert aesthetics.title == "Scale"
        assert aesthetics.font_size == 8

    def test_custom_values(self):
        """Test custom values."""
        aesthetics = ScaleBarAesthetics(
            style="std",
            color="blue",
            position="center",
            scale_distance=500,
            title="Custom Scale",
        )

        assert aesthetics.style == "std"
        assert aesthetics.color == "blue"
        assert aesthetics.position == "center"
        assert aesthetics.scale_distance == 500
        assert aesthetics.title == "Custom Scale"

    def test_invalid_position(self):
        """Test invalid position value."""
        with pytest.raises(ValidationError):
            ScaleBarAesthetics(position="invalid_position")


class TestGenesAesthetics:
    """Test cases for the GenesAesthetics class."""

    def test_default_values(self):
        """Test default values."""
        aesthetics = GenesAesthetics()

        assert aesthetics.style == "std"
        assert aesthetics.color == "#2c3e50"  # Genome browser dark navy
        assert aesthetics.fill is True
        assert aesthetics.alpha == 1.0
        assert aesthetics.display == "collapsed"
        assert aesthetics.minimum_gene_length == 0
        assert aesthetics.max_number_of_rows == 4
        assert aesthetics.interval_height == 0.1  # Further reduced

    def test_custom_values(self):
        """Test custom values."""
        aesthetics = GenesAesthetics(
            style="std",
            color="green",
            fill=False,
            alpha=0.7,
            display="expanded",
            minimum_gene_length=500,
            max_number_of_rows=10,
            interval_height=0.3,
        )

        assert aesthetics.style == "std"
        assert aesthetics.color == "green"
        assert aesthetics.fill is False
        assert aesthetics.alpha == 0.7
        assert aesthetics.display == "expanded"
        assert aesthetics.minimum_gene_length == 500
        assert aesthetics.max_number_of_rows == 10
        assert aesthetics.interval_height == 0.3

    def test_invalid_display(self):
        """Test invalid display value."""
        with pytest.raises(ValidationError):
            GenesAesthetics(display="invalid_display")
