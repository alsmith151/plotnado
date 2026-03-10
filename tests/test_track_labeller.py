"""
Tests for the TrackLabeller class.
"""

import pytest
from unittest.mock import patch, call, MagicMock

from plotnado.tracks import (
    TrackLabeller,
    GenomicRegion,
    clean_axis
)


class TestTrackLabeller:
    """Test cases for TrackLabeller class."""

    def test_init(self):
        """Test initializing a TrackLabeller instance."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0
        )
        
        assert labeller.gr == gr
        assert labeller.y_min == 0.0
        assert labeller.y_max == 1.0
        assert labeller.title == "Test Title"
        assert labeller.scale_min == 0.0
        assert labeller.scale_max == 1.0
        assert labeller.title_location == "left"
        assert labeller.scale_location == "right"

    def test_y_delta_property(self):
        """Test the y_delta property."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0
        )
        
        assert labeller.y_delta == 1.0

    def test_format_scale(self):
        """Test the _format_scale method."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0,
            scale_precision=2
        )
        
        # Test integer values
        assert labeller._format_scale(5) == "5"
        assert labeller._format_scale(10.0) == "10"
        
        # Test float values
        assert labeller._format_scale(0.123) == "0.12"
        assert labeller._format_scale(0.126) == "0.13"  # Rounds

    def test_plot_title_left(self, mock_ax):
        """Test plotting title at left position."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0,
            title_location="left"
        )
        
        labeller._plot_title(mock_ax, gr)
        
        # Validate
        mock_ax.text.assert_called_once()
        args, kwargs = mock_ax.text.call_args
        assert args[0] == 1010  # gr.start + (0.01 * gr.length)
        assert args[1] == 0.8  # labeller.y_delta * labeller.title_height
        assert args[2] == "Test Title"
        assert kwargs['horizontalalignment'] == "left"
        assert kwargs['verticalalignment'] == "top"
        
    def test_plot_title_right(self, mock_ax):
        """Test plotting title at right position."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0,
            title_location="right"
        )
        
        labeller._plot_title(mock_ax, gr)
        
        # Validate
        mock_ax.text.assert_called_once()
        args, kwargs = mock_ax.text.call_args
        assert args[0] == 1990  # gr.end - (0.01 * gr.length)
        assert args[1] == 0.8  # labeller.y_delta * labeller.title_height
        assert args[2] == "Test Title"
        assert kwargs['horizontalalignment'] == "right"
        assert kwargs['verticalalignment'] == "top"

    def test_plot_scale_left(self, mock_ax):
        """Test plotting scale at left position."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.5,
            y_max=1.5,
            title="Test Title",
            scale_min=0.5,
            scale_max=1.5,
            scale_location="left"
        )
        
        labeller._plot_scale(mock_ax, gr)
        
        # Validate
        mock_ax.text.assert_called_once()
        args, kwargs = mock_ax.text.call_args
        assert args[0] == 1010  # gr.start + (0.01 * gr.length)
        assert args[1] == 0.8  # labeller.y_delta * labeller.scale_height
        assert args[2] == "[ 0.50 - 1.50 ]"
        assert kwargs['horizontalalignment'] == "left"
        assert kwargs['verticalalignment'] == "top"

    def test_plot_scale_right(self, mock_ax):
        """Test plotting scale at right position."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.5,
            y_max=1.5,
            title="Test Title",
            scale_min=0.5,
            scale_max=1.5,
            scale_location="right"
        )
        
        labeller._plot_scale(mock_ax, gr)
        
        # Validate
        mock_ax.text.assert_called_once()
        args, kwargs = mock_ax.text.call_args
        assert args[0] == 1990  # gr.end - (0.01 * gr.length)
        assert args[1] == 0.8  # labeller.y_delta * labeller.scale_height
        assert args[2] == "[ 0.50 - 1.50 ]"
        assert kwargs['horizontalalignment'] == "right"
        assert kwargs['verticalalignment'] == "top"

    def test_label_box_enabled_default(self, mock_ax):
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")

        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            plot_title=True,
            plot_scale=False,
        )

        labeller.plot(mock_ax, gr)
        _, kwargs = mock_ax.text.call_args
        assert kwargs["bbox"]["facecolor"] == "white"
        assert kwargs["bbox"]["alpha"] == 0.9

    def test_label_box_disabled(self, mock_ax):
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")

        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            plot_title=True,
            plot_scale=False,
            label_box_enabled=False,
        )

        labeller.plot(mock_ax, gr)
        _, kwargs = mock_ax.text.call_args
        assert kwargs["bbox"] is None

    def test_label_box_alpha_custom(self, mock_ax):
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")

        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            plot_title=True,
            plot_scale=False,
            label_box_alpha=0.5,
        )

        labeller.plot(mock_ax, gr)
        _, kwargs = mock_ax.text.call_args
        assert kwargs["bbox"]["alpha"] == 0.5

    @patch.object(TrackLabeller, '_plot_title')
    @patch.object(TrackLabeller, '_plot_scale')
    @patch('plotnado.tracks.clean_axis')
    def test_plot_all_enabled(self, mock_clean_axis, mock_plot_scale, mock_plot_title, mock_ax):
        """Test plotting with all elements enabled."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0,
            plot_title=True,
            plot_scale=True
        )
        
        labeller.plot(mock_ax, gr)
        
        # Validate
        mock_plot_title.assert_called_once_with(mock_ax, gr)
        mock_plot_scale.assert_called_once_with(mock_ax, gr)
        mock_clean_axis.assert_called_once_with(mock_ax)

    @patch('plotnado.tracks.clean_axis')
    def test_plot_forces_scale_opposite_title_left(self, mock_clean_axis, mock_ax):
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")

        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            title_location="left",
            scale_location="left",  # intentionally same as title
            plot_title=True,
            plot_scale=True,
        )

        labeller.plot(mock_ax, gr)

        # 1st text call is title, 2nd is scale
        scale_call = mock_ax.text.call_args_list[1]
        args, kwargs = scale_call
        assert args[0] == 1990
        assert kwargs["horizontalalignment"] == "right"
        mock_clean_axis.assert_called_once_with(mock_ax)

    @patch('plotnado.tracks.clean_axis')
    def test_plot_forces_scale_opposite_title_right(self, mock_clean_axis, mock_ax):
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")

        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            title_location="right",
            scale_location="right",  # intentionally same as title
            plot_title=True,
            plot_scale=True,
        )

        labeller.plot(mock_ax, gr)

        scale_call = mock_ax.text.call_args_list[1]
        args, kwargs = scale_call
        assert args[0] == 1010
        assert kwargs["horizontalalignment"] == "left"
        mock_clean_axis.assert_called_once_with(mock_ax)

    @patch('plotnado.tracks.clean_axis')
    def test_plot_keeps_scale_location_when_no_title(self, mock_clean_axis, mock_ax):
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")

        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="",
            scale_location="left",
            plot_title=True,
            plot_scale=True,
        )

        labeller.plot(mock_ax, gr)

        args, kwargs = mock_ax.text.call_args
        assert args[0] == 1010
        assert kwargs["horizontalalignment"] == "left"
        mock_clean_axis.assert_called_once_with(mock_ax)

    @patch.object(TrackLabeller, '_plot_title')
    @patch.object(TrackLabeller, '_plot_scale')
    @patch('plotnado.tracks.clean_axis')
    def test_plot_title_only(self, mock_clean_axis, mock_plot_scale, mock_plot_title, mock_ax):
        """Test plotting with only title enabled."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0,
            plot_title=True,
            plot_scale=False
        )
        
        labeller.plot(mock_ax, gr)
        
        # Validate
        mock_plot_title.assert_called_once_with(mock_ax, gr)
        mock_plot_scale.assert_not_called()
        mock_clean_axis.assert_called_once_with(mock_ax)

    @patch.object(TrackLabeller, '_plot_title')
    @patch.object(TrackLabeller, '_plot_scale')
    @patch('plotnado.tracks.clean_axis')
    def test_plot_scale_only(self, mock_clean_axis, mock_plot_scale, mock_plot_title, mock_ax):
        """Test plotting with only scale enabled."""
        gr = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        labeller = TrackLabeller(
            gr=gr,
            y_min=0.0,
            y_max=1.0,
            title="Test Title",
            scale_min=0.0,
            scale_max=1.0,
            plot_title=False,
            plot_scale=True
        )
        
        labeller.plot(mock_ax, gr)
        
        # Validate
        mock_plot_title.assert_not_called()
        mock_plot_scale.assert_called_once_with(mock_ax, gr)
        mock_clean_axis.assert_called_once_with(mock_ax)


class TestCleanAxis:
    """Test cases for clean_axis function."""

    def test_clean_axis(self, mock_ax):
        """Test clean_axis function."""
        clean_axis(mock_ax)
        
        # Validate
        mock_ax.set_xticks.assert_called_once_with([])
        mock_ax.set_yticks.assert_called_once_with([])
        mock_ax.spines["top"].set_visible.assert_called_once_with(False)
        mock_ax.spines["right"].set_visible.assert_called_once_with(False)
        mock_ax.spines["left"].set_visible.assert_called_once_with(False)
        mock_ax.spines["bottom"].set_visible.assert_called_once_with(False)
        mock_ax.xaxis.set_major_locator.assert_called_once()
        mock_ax.yaxis.set_major_locator.assert_called_once()
