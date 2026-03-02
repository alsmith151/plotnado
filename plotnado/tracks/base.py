"""
Base track classes.
"""

from abc import ABC
from typing import Any, Literal

import matplotlib.axes
from pydantic import BaseModel, ConfigDict
from .region import GenomicRegion


class Track(BaseModel, ABC):
    """
    Abstract base class for all track types.

    Attributes:
        title: Optional title for the track
        height: Relative height of the track (default 1.0)
    """

    title: str | None = None
    data: Any | None = None
    height: float = 1.0
    autoscale_group: str | None = None
    label_on_track: bool = False
    data_range_style: Literal["text", "colorbar", "none"] = "text"
    label_box_enabled: bool = True
    label_box_alpha: float = 0.9

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def fetch_data(self, gr: GenomicRegion) -> Any:
        """Fetch data for the given genomic region."""
        raise NotImplementedError

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the track on the given axes for the given region."""
        raise NotImplementedError


class TrackLabeller(BaseModel):
    """
    Handles track labelling (title and scale display).
    """

    gr: GenomicRegion
    y_min: float
    y_max: float

    plot_title: bool = True
    plot_scale: bool = True
    label_on_track: bool = False
    data_range_style: Literal["text", "colorbar", "none"] = "text"
    label_box_enabled: bool = True
    label_box_alpha: float = 0.9

    title: str = ""
    title_size: int = 10
    title_color: str = "#333333"  # Darker gray
    title_font: str = "DejaVu Sans"
    title_weight: Literal["normal", "bold"] = "bold"
    title_location: Literal["left", "right"] = "left"
    title_height: float = 0.8

    scale_min: float = 0
    scale_max: float = 1
    scale_precision: int = 2
    scale_size: int = 9
    scale_color: str = "#666666"  # Gray
    scale_font: str = "DejaVu Sans"
    scale_weight: Literal["normal", "bold"] = "normal"
    scale_location: Literal["left", "right"] = "right"
    scale_height: float = 0.8

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @staticmethod
    def _default_text_bbox(alpha: float) -> dict:
        return {
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": alpha,
            "boxstyle": "round,pad=0.2",
        }

    def _text_bbox(self) -> dict | None:
        if not self.label_box_enabled:
            return None
        return self._default_text_bbox(self.label_box_alpha)

    @property
    def y_delta(self) -> float:
        return self.y_max - self.y_min

    def _plot_title(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        if self.label_on_track:
            x_pos = gr.start + (0.02 * gr.length)
            y_pos = self.y_min + (self.y_delta * self.title_height)
            h_align = "left"
        else:
            x_pos = (
                gr.start + (0.01 * gr.length)
                if self.title_location == "left"
                else gr.end - (0.01 * gr.length)
            )
            y_pos = self.y_delta * self.title_height
            h_align = "left" if self.title_location == "left" else "right"

        ax.text(
            x_pos,
            y_pos,
            self.title,
            horizontalalignment=h_align,
            verticalalignment="top",
            bbox=self._text_bbox(),
            fontdict={
                "size": self.title_size,
                "color": self.title_color,
                "fontname": self.title_font,
                "weight": self.title_weight,
            },
        )

    def _format_scale(self, value: float) -> str:
        """Format scale value for display."""
        # Use integer if it's a whole number
        if value % 1 == 0:
            return str(int(value))
        return f"{value:.{self.scale_precision}f}"

    def _plot_scale(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        self._plot_scale_at(ax, gr, self.scale_location)

    def _plot_scale_at(
        self,
        ax: matplotlib.axes.Axes,
        gr: GenomicRegion,
        location: Literal["left", "right"],
    ) -> None:
        y_min = self._format_scale(self.y_min)
        y_max = self._format_scale(self.y_max)

        x_pos = (
            gr.end - (0.01 * gr.length)
            if location == "right"
            else gr.start + (0.01 * gr.length)
        )
        h_align = "right" if location == "right" else "left"

        ax.text(
            x_pos,
            self.y_delta * self.scale_height,
            f"[ {y_min} - {y_max} ]",
            horizontalalignment=h_align,
            verticalalignment="top",
            bbox=self._text_bbox(),
            fontdict={
                "size": self.scale_size,
                "color": self.scale_color,
                "fontname": self.scale_font,
                "weight": self.scale_weight,
            },
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> "TrackLabeller":
        should_plot_scale = self.plot_scale and self.data_range_style == "text"
        should_plot_title = self.plot_title and bool(self.title)

        title_side: Literal["left", "right"]
        if self.label_on_track:
            title_side = "left"
        else:
            title_side = self.title_location

        if should_plot_title:
            self._plot_title(ax, gr)
        if should_plot_scale:
            original_scale_location = self.scale_location
            if should_plot_title:
                self.scale_location = "right" if title_side == "left" else "left"
            try:
                self._plot_scale(ax, gr)
            finally:
                self.scale_location = original_scale_location

        import plotnado.tracks as pg_tracks

        pg_tracks.clean_axis(ax)
        return self
