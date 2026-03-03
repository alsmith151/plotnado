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

    def test_flattened_track_kwargs_sync_to_aesthetics(self):
        """Track constructor should accept flattened aesthetics kwargs."""
        track = BigWigTrack(color="red", alpha=0.35, style="line")

        assert track.color == "red"
        assert track.alpha == 0.35
        assert track.style == "line"
        assert track.aesthetics.color == "red"
        assert track.aesthetics.alpha == 0.35
        assert track.aesthetics.style == "line"

    def test_nested_aesthetics_still_works(self):
        """Legacy nested aesthetics input remains supported."""
        track = BigWigTrack(aesthetics=BigwigAesthetics(color="purple", alpha=0.5))

        assert track.color == "purple"
        assert track.alpha == 0.5
        assert track.aesthetics.color == "purple"
        assert track.aesthetics.alpha == 0.5

    def test_explicit_flattened_kwargs_override_nested(self):
        """Explicit flattened kwargs should win over nested aesthetics."""
        track = BigWigTrack(
            color="orange",
            aesthetics=BigwigAesthetics(color="purple", alpha=0.5),
        )

        assert track.color == "orange"
        assert track.aesthetics.color == "orange"

    def test_signature_exposes_flattened_aesthetic_fields(self):
        """Constructor signature should expose flattened kwargs for notebook discoverability."""
        signature = inspect.signature(BigWigTrack)

        assert "color" in signature.parameters
        assert "alpha" in signature.parameters
        assert "style" in signature.parameters

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

    def test_default_values(self):
        """Test default values."""
        aesthetics = ScaleBarAesthetics()

        assert aesthetics.style == "std"
        assert aesthetics.color == "#333333"  # Dark gray
        assert aesthetics.position == "left"
        assert aesthetics.scale_distance is None
        assert aesthetics.title == "Scale"
        assert aesthetics.font_size == 8
        assert aesthetics.bar_linewidth == 3.0
        assert aesthetics.tick_linewidth == 2.0
        assert aesthetics.tick_height == 0.1
        assert aesthetics.label_offset == 0.25

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
        assert aesthetics.color == "black"
        assert aesthetics.fill is True
        assert aesthetics.alpha == 1.0
        assert aesthetics.display == "collapsed"
        assert aesthetics.minimum_gene_length == 0
        assert aesthetics.max_number_of_rows == 4
        assert aesthetics.interval_height == 0.1  # Further reduced
        assert aesthetics.exon_linewidth == 0.8
        assert aesthetics.exon_edge_color == "black"
        assert aesthetics.exon_color == "black"
        assert aesthetics.intron_linewidth == 0.8
        assert aesthetics.intron_color == "black"
        assert aesthetics.gene_label_font_size == 8
        assert aesthetics.gene_label_style == "italic"

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
    def test_default_values(self):
        aesthetics = GenomicAxisAesthetics()

        assert aesthetics.color == "#666666"
        assert aesthetics.font_size == 9
        assert aesthetics.num_ticks == 5
        assert aesthetics.show_chromosome is True
        assert aesthetics.tick_height == 0.15
        assert aesthetics.axis_linewidth == 1.5
        assert aesthetics.tick_color == "#333333"
        assert aesthetics.tick_linewidth == 1.2
        assert aesthetics.chromosome_fontweight == "bold"

    def test_genomic_axis_flattens_aesthetics(self):
        track = GenomicAxis(axis_linewidth=2.0, tick_color="purple")

        assert track.axis_linewidth == 2.0
        assert track.tick_color == "purple"
        assert track.aesthetics.axis_linewidth == 2.0
        assert track.aesthetics.tick_color == "purple"


class TestAdditionalAesthetics:
    def test_bed_defaults(self):
        aesthetics = BedAesthetics()
        assert aesthetics.rect_linewidth == 1.0

    def test_highlight_defaults(self):
        aesthetics = HighlightsAesthetics()
        assert aesthetics.linewidth == 1.0

    def test_links_defaults(self):
        aesthetics = LinksAesthetics()
        assert aesthetics.y_baseline == 0.1
