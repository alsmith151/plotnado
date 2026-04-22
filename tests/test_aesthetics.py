"""
Tests for the aesthetics classes.
"""

import inspect

import pytest
from pydantic import ValidationError

from plotnado.tracks import (
    BedAesthetics,
    BigwigAesthetics,
    BigWigTrack,
    GenomicAxis,
    GenomicAxisAesthetics,
    HighlightsAesthetics,
    LabelConfig,
    LinksAesthetics,
    list_options,
    ScaleBarAesthetics,
    GenesAesthetics,
)


class TestBigwigAesthetics:
    """Test cases for the BigwigAesthetics class."""

    def test_instantiates_with_defaults(self):
        aesthetics = BigwigAesthetics()
        assert aesthetics.style is not None
        assert aesthetics.color is not None

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
        )

        assert aesthetics.style == "scatter"
        assert aesthetics.color == "red"
        assert aesthetics.fill is False
        assert aesthetics.alpha == 0.5
        assert aesthetics.scatter_point_size == 2.0
        assert aesthetics.linewidth == 2.0
        assert aesthetics.min_value == 0.0
        assert aesthetics.max_value == 1.0

    def test_invalid_style(self):
        """Test invalid style value."""
        with pytest.raises(ValidationError):
            BigwigAesthetics(style="invalid_style")

    def test_track_uses_nested_aesthetics(self):
        """Track constructor accepts nested aesthetics only."""
        track = BigWigTrack(aesthetics=BigwigAesthetics(color="purple", alpha=0.5))

        assert track.color == "purple"
        assert track.alpha == 0.5
        assert track.aesthetics.color == "purple"
        assert track.aesthetics.alpha == 0.5

    def test_track_signature_does_not_expose_flattened_aesthetics(self):
        """Constructor signature should expose nested aesthetics model only."""
        signature = inspect.signature(BigWigTrack)

        assert "aesthetics" in signature.parameters
        assert "color" not in signature.parameters
        assert "alpha" not in signature.parameters
        assert "style" not in signature.parameters

    def test_track_options_reports_track_and_aesthetic_sections(self):
        """Programmatic option discovery should split track vs aesthetics fields."""
        options = BigWigTrack.options()

        assert "track" in options
        assert "aesthetics" in options
        assert "label" in options
        assert "title" in options["track"]
        assert "height" in options["track"]
        assert "color" in options["aesthetics"]
        assert "alpha" in options["aesthetics"]
        assert "style" in options["aesthetics"]
        assert options["aesthetics"]["style"]["choices"] == [
            "std",
            "fill",
            "line",
            "scatter",
            "heatmap",
            "fragment",
        ]
        assert "plot_title" in options["label"]

    def test_list_options_helper_matches_track_options(self):
        """Module helper should mirror Track.options for notebook workflows."""
        assert list_options(BigWigTrack) == BigWigTrack.options()

    def test_options_markdown_contains_sections(self):
        markdown = BigWigTrack.options_markdown()

        assert "## BigWigTrack options" in markdown
        assert "### Track fields" in markdown
        assert "### Aesthetics fields" in markdown
        assert "### Label fields" in markdown
        assert "| color |" in markdown

    def test_track_label_config(self):
        track = BigWigTrack(label=LabelConfig(plot_title=False, plot_scale=False))
        assert track.label.plot_title is False
        assert track.label.plot_scale is False


class TestScaleBarAesthetics:
    """Test cases for the ScaleBarAesthetics class."""

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


class TestGenomicAxisAesthetics:
    def test_genomic_axis_uses_nested_aesthetics(self):
        track = GenomicAxis(
            aesthetics=GenomicAxisAesthetics(axis_linewidth=2.0, tick_color="purple")
        )

        assert track.axis_linewidth == 2.0
        assert track.tick_color == "purple"
        assert track.aesthetics.axis_linewidth == 2.0
        assert track.aesthetics.tick_color == "purple"


class TestAdditionalAesthetics:
    def test_bed_instantiates(self):
        BedAesthetics()

    def test_highlight_instantiates(self):
        HighlightsAesthetics()

    def test_links_instantiates(self):
        LinksAesthetics()

    def test_narrowpeak_instantiates(self):
        from plotnado.tracks import NarrowPeakAesthetics

        NarrowPeakAesthetics()
