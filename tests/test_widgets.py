"""Tests for notebook track visibility widgets."""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")

from plotnado import GenomicFigure, TrackVisibilityWidget


def test_track_visibility_widget_toggles_visible_track(bedgraph_df) -> None:
    fig = GenomicFigure()
    fig.bigwig(bedgraph_df, title="Signal")

    widget = fig.track_visibility_widget("chr1:100-600")

    assert isinstance(widget, TrackVisibilityWidget)
    signal_control = next(control for control in widget.controls if control.label == "Signal")
    assert signal_control.checkbox.value is True

    signal_control.checkbox.value = False
    assert fig["Signal"].height == 0.0

    signal_control.checkbox.value = True
    assert fig["Signal"].height == 1.0


def test_track_visibility_widget_restores_hidden_track_and_skips_overlays(
    bedgraph_df,
    bed_df,
) -> None:
    fig = GenomicFigure()
    fig.bigwig(bedgraph_df, title="Signal")
    fig.bed(bed_df, title="Hidden bed", height=0.0)
    fig.vline(250)

    widget = fig.track_visibility_widget("chr1:100-600")

    labels = [control.label for control in widget.controls]
    assert "Hidden bed" in labels
    assert "VLineTrack" not in labels

    hidden_control = next(control for control in widget.controls if control.label == "Hidden bed")
    assert hidden_control.checkbox.value is False

    hidden_control.checkbox.value = True
    assert fig["Hidden bed"].height == 1.0