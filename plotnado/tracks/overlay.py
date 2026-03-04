"""
BigWig overlay track for displaying multiple signals on the same axis.
"""

from pathlib import Path

import matplotlib.axes
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from .region import GenomicRegion
from .base import LabelConfig, Track, TrackLabeller
from .utils import clean_axis
from .bigwig import BigWigTrack, BigwigAesthetics
from .scaling import Autoscaler


class BigwigOverlayAesthetics(BaseModel):
    """
    Aesthetics for BigWig overlay tracks.
    """

    colors: list[str] | None = Field(
        default=None,
        description="Optional per-track color overrides for overlaid subtracks.",
    )
    alpha: float = Field(default=0.5, description="Opacity applied to overlaid signal traces.")
    show_labels: bool = Field(default=True, description="Render labels for component tracks in overlay contexts.")
    min_value: float | None = Field(
        default=None,
        description="Optional shared minimum y-value for all overlaid tracks.",
    )
    max_value: float | None = Field(
        default=None,
        description="Optional shared maximum y-value for all overlaid tracks.",
    )


class BigwigOverlay(Track):
    """
    Track for overlaying multiple BigWig signals.

    Attributes:
        tracks: List of BigWigTrack instances or paths
        aesthetics: Visual configuration
    """

    tracks: list[BigWigTrack | Path | str] = Field(
        description="List of BigWigTrack instances or BigWig file paths to overlay.",
    )
    aesthetics: BigwigOverlayAesthetics = Field(
        default_factory=BigwigOverlayAesthetics,
        description="Overlay-level visual styling and shared y-range settings.",
    )
    height: float = Field(default=2.0, description="Relative panel height for this track.")

    _track_instances: list[BigWigTrack] = []

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __init__(self, **data):
        super().__init__(**data)
        self._track_instances = []
        default_colors = ["steelblue", "darkorange", "forestgreen", "crimson", "purple"]

        for i, t in enumerate(self.tracks):
            if isinstance(t, BigWigTrack):
                inst = t
            else:
                color = (
                    self.colors[i]
                    if self.colors and i < len(self.colors)
                    else default_colors[i % len(default_colors)]
                )
                inst = BigWigTrack(
                    data=t,
                    title=Path(t).stem if isinstance(t, (Path, str)) else f"Track {i}",
                    aesthetics=BigwigAesthetics(
                        color=color,
                        alpha=self.alpha,
                    ),
                    label=LabelConfig(plot_title=False, plot_scale=False),
                )
            self._track_instances.append(inst)

    def fetch_data(self, gr: GenomicRegion) -> list[pd.DataFrame]:
        """Fetch data from all subtracks."""
        return [t.fetch_data(gr) for t in self._track_instances]

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot overlaid tracks."""
        # Calculate global scale if auto
        scaler = Autoscaler(tracks=self._track_instances, gr=gr)
        y_min = (
            self.min_value
            if self.min_value is not None
            else scaler.min_value
        )
        y_max = (
            self.max_value
            if self.max_value is not None
            else scaler.max_value
        )

        # Ensure some range
        if y_min == y_max:
            y_max = y_min + 1

        for track in self._track_instances:
            # Temporarily override aesthetics for plotting
            orig_min = track.min_value
            orig_max = track.max_value
            orig_scale = track.label.plot_scale
            orig_title = track.label.plot_title

            track.min_value = y_min
            track.max_value = y_max
            track.label.plot_scale = False
            track.label.plot_title = False

            track.plot(ax, gr)

            # Restore
            track.min_value = orig_min
            track.max_value = orig_max
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
            )
            labeller.plot(ax, gr)
        else:
            clean_axis(ax)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
