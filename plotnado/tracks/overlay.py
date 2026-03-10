"""Generic overlay track for rendering multiple tracks on the same axis."""

from pathlib import Path

import matplotlib.axes
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr

from .region import GenomicRegion
from .base import LabelConfig, Track, TrackLabeller
from .utils import clean_axis
from .bigwig import BigWigTrack, BigwigAesthetics


class OverlayTrackAesthetics(BaseModel):
    """
    Aesthetics for overlay tracks.
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
                inst = BigWigTrack(
                    data=t,
                    title=Path(t).stem if isinstance(t, (Path, str)) else f"Track {i}",
                    aesthetics=BigwigAesthetics(
                        color=color,
                        alpha=self.alpha,
                    ),
                    label=LabelConfig(plot_title=False, plot_scale=False),
                )

            if self.colors and i < len(self.colors) and inst.has_aesthetic("color"):
                inst.color = self.colors[i]
            if inst.has_aesthetic("alpha"):
                inst.alpha = self.alpha

            self._track_instances.append(inst)

    def fetch_data(self, gr: GenomicRegion) -> list[pd.DataFrame]:
        """Fetch data from all subtracks."""
        return [t.fetch_data(gr) for t in self._track_instances]

    def _shared_limits(self, gr: GenomicRegion) -> tuple[float, float]:
        values: list[np.ndarray] = []

        for track in self._track_instances:
            data = track.fetch_data(gr)
            if isinstance(data, pd.DataFrame):
                if "value" in data.columns:
                    values.append(data["value"].to_numpy(dtype=float))
                elif not data.empty:
                    candidate = pd.to_numeric(data.iloc[:, -1], errors="coerce").to_numpy(dtype=float)
                    values.append(candidate)
            elif isinstance(data, np.ndarray):
                values.append(data.astype(float).ravel())

        if not values:
            return 0.0, 1.0

        merged = np.concatenate(values)
        merged = merged[~np.isnan(merged)]
        if merged.size == 0:
            return 0.0, 1.0

        y_min = self.min_value if self.min_value is not None else float(min(0.0, np.min(merged)))
        y_max = self.max_value if self.max_value is not None else float(np.max(merged))
        if y_min == y_max:
            y_max = y_min + 1.0
        return y_min, y_max

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
