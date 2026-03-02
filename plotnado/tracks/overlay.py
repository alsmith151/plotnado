"""
BigWig overlay track for displaying multiple signals on the same axis.
"""

from pathlib import Path

import matplotlib.axes
import pandas as pd
from pydantic import BaseModel, ConfigDict

from .region import GenomicRegion
from .base import Track, TrackLabeller
from .utils import clean_axis
from .bigwig import BigWigTrack, BigwigAesthetics
from .scaling import Autoscaler


class BigwigOverlayAesthetics(BaseModel):
    """
    Aesthetics for BigWig overlay tracks.
    """

    colors: list[str] | None = None
    alpha: float = 0.5
    show_labels: bool = True
    min_value: float | None = None
    max_value: float | None = None


class BigwigOverlay(Track):
    """
    Track for overlaying multiple BigWig signals.

    Attributes:
        tracks: List of BigWigTrack instances or paths
        aesthetics: Visual configuration
    """

    tracks: list[BigWigTrack | Path | str]
    aesthetics: BigwigOverlayAesthetics = BigwigOverlayAesthetics()
    height: float = 2.0

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
                    self.aesthetics.colors[i]
                    if self.aesthetics.colors and i < len(self.aesthetics.colors)
                    else default_colors[i % len(default_colors)]
                )
                inst = BigWigTrack(
                    data=t,
                    title=Path(t).stem if isinstance(t, (Path, str)) else f"Track {i}",
                    aesthetics=BigwigAesthetics(
                        color=color,
                        alpha=self.aesthetics.alpha,
                        plot_title=False,
                        plot_scale=False,
                    ),
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
            self.aesthetics.min_value
            if self.aesthetics.min_value is not None
            else scaler.min_value
        )
        y_max = (
            self.aesthetics.max_value
            if self.aesthetics.max_value is not None
            else scaler.max_value
        )

        # Ensure some range
        if y_min == y_max:
            y_max = y_min + 1

        for track in self._track_instances:
            # Temporarily override aesthetics for plotting
            orig_min = track.aesthetics.min_value
            orig_max = track.aesthetics.max_value
            orig_scale = track.aesthetics.plot_scale
            orig_title = track.aesthetics.plot_title

            track.aesthetics.min_value = y_min
            track.aesthetics.max_value = y_max
            track.aesthetics.plot_scale = False
            track.aesthetics.plot_title = False

            track.plot(ax, gr)

            # Restore
            track.aesthetics.min_value = orig_min
            track.aesthetics.max_value = orig_max
            track.aesthetics.plot_scale = orig_scale
            track.aesthetics.plot_title = orig_title

        # Add global labeller
        if self.title:
            labeller = TrackLabeller(
                gr=gr,
                y_min=y_min,
                y_max=y_max,
                title=self.title,
                plot_title=True,
                plot_scale=True,
                scale_min=y_min,
                scale_max=y_max,
                label_box_enabled=self.label_box_enabled,
                label_box_alpha=self.label_box_alpha,
            )
            labeller.plot(ax, gr)
        else:
            clean_axis(ax)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
