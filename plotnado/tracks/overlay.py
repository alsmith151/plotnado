"""Generic overlay track for rendering multiple tracks on the same axis."""

from pathlib import Path

import matplotlib.axes
import numpy as np
import pandas as pd
from pydantic import ConfigDict, Field, PrivateAttr

from .region import GenomicRegion
from .base import LabelConfig, Track, TrackLabeller
from .scaling import calculate_data_limits
from .utils import clean_axis
from .bigwig import BigWigTrack, BigwigAesthetics
from .enums import PlotStyle, TrackType
from .aesthetics import BaseMultiColorAesthetics
from .registry import registry


class OverlayTrackAesthetics(BaseMultiColorAesthetics):
    """
    Aesthetics for overlay tracks.

    Inherits colors and alpha from BaseMultiColorAesthetics.
    """

    show_labels: bool = Field(default=True, description="Render labels for component tracks in overlay contexts.")
    style: PlotStyle | None = Field(
        default=None,
        description="Optional BigWig render style applied to component tracks that support it.",
    )
    min_value: float | None = Field(
        default=None,
        description="Optional shared minimum y-value for all overlaid tracks.",
    )
    max_value: float | None = Field(
        default=None,
        description="Optional shared maximum y-value for all overlaid tracks.",
    )


@registry.register(TrackType.OVERLAY, aliases=["bigwig_overlay"])
class OverlayTrack(Track):
    """
    Track for overlaying multiple tracks on the same axis.

    Attributes:
        tracks: List of Track instances or file-like BigWig inputs
        aesthetics: Visual configuration
    """

    tracks: list[Track | pd.DataFrame | Path | str] = Field(
        description="List of Track instances, in-memory signal DataFrames, or BigWig-like file paths to overlay.",
    )
    aesthetics: OverlayTrackAesthetics = Field(
        default_factory=OverlayTrackAesthetics,
        description="Overlay-level visual styling and shared y-range settings.",
    )
    height: float = Field(default=2.0, description="Relative panel height for this track.")

    _track_instances: list[Track] = PrivateAttr(default_factory=list)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __init__(self, **data):
        super().__init__(**data)
        self._track_instances = []
        default_colors = ["steelblue", "darkorange", "forestgreen", "crimson", "purple"]

        for i, t in enumerate(self.tracks):
            if isinstance(t, Track):
                inst = t
            else:
                color = (
                    self.colors[i]
                    if self.colors and i < len(self.colors)
                    else default_colors[i % len(default_colors)]
                )
                aesthetics_kwargs = {
                    "color": color,
                    "alpha": self.alpha,
                }
                if self.style is not None:
                    aesthetics_kwargs["style"] = self.style
                inst = BigWigTrack(
                    data=t,
                    title=Path(t).stem if isinstance(t, (Path, str)) else f"Track {i}",
                    aesthetics=BigwigAesthetics(**aesthetics_kwargs),
                    label=LabelConfig(plot_title=False, plot_scale=False),
                )

            if self.colors and i < len(self.colors) and inst.has_aesthetic("color"):
                inst.color = self.colors[i]
            if inst.has_aesthetic("alpha"):
                inst.alpha = self.alpha
            if self.style is not None and inst.has_aesthetic("style"):
                inst.style = self.style

            self._track_instances.append(inst)

    def fetch_data(self, gr: GenomicRegion) -> list[pd.DataFrame]:
        """Fetch data from all subtracks."""
        return [t.fetch_data(gr) for t in self._track_instances]

    def _shared_limits(self, gr: GenomicRegion) -> tuple[float, float]:
        return calculate_data_limits(
            self.fetch_data(gr),
            min_value=self.min_value,
            max_value=self.max_value,
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot overlaid tracks."""
        y_min, y_max = self._shared_limits(gr)

        for track in self._track_instances:
            # Temporarily override aesthetics for plotting
            orig_min = getattr(track, "min_value", None)
            orig_max = getattr(track, "max_value", None)
            has_label = hasattr(track, "label") and track.label is not None
            if has_label:
                orig_scale = track.label.plot_scale
                orig_title = track.label.plot_title

            if track.has_aesthetic("min_value"):
                track.min_value = y_min
            if track.has_aesthetic("max_value"):
                track.max_value = y_max
            if has_label and not self.show_labels:
                track.label.plot_scale = False
                track.label.plot_title = False

            track.plot(ax, gr)

            # Restore
            if track.has_aesthetic("min_value"):
                track.min_value = orig_min
            if track.has_aesthetic("max_value"):
                track.max_value = orig_max
            if has_label and not self.show_labels:
                track.label.plot_scale = orig_scale
                track.label.plot_title = orig_title

        # Add global labeller
        if self.title:
            labeller = TrackLabeller.from_config(
                self.label,
                gr,
                y_min,
                y_max,
                title=self.title,
                title_color=self._track_instances[0].color if self._track_instances else None,
            )
            labeller.plot(ax, gr)
        else:
            clean_axis(ax)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)


class BigwigOverlay(OverlayTrack):
    """Backward-compatible alias for OverlayTrack."""


BigwigOverlayAesthetics = OverlayTrackAesthetics
