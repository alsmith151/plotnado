"""
Base track classes.
"""

from abc import ABC
from typing import Any, Literal, Optional

import matplotlib.axes
from pydantic import BaseModel
from .region import GenomicRegion


class Track(BaseModel, ABC):
    """
    Abstract base class for all track types.

    Attributes:
        title: Optional title for the track
        height: Relative height of the track (default 1.0)
    """

    title: Optional[str] = None
    data: Optional[Any] = None
    height: float = 1.0

    class Config:
        arbitrary_types_allowed = True

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

    title: str = ""
    title_size: int = 8
    title_color: str = "black"
    title_font: str = "Arial"
    title_weight: Literal["normal", "bold"] = "bold"
    title_location: Literal["left", "right"] = "left"
    title_height: float = 0.8

    scale_min: float = 0
    scale_max: float = 1
    scale_precision: int = 2
    scale_size: int = 6
    scale_color: str = "black"
    scale_font: str = "Arial"
    scale_weight: Literal["normal", "bold"] = "bold"
    scale_location: Literal["left", "right"] = "right"
    scale_height: float = 0.8

    class Config:
        arbitrary_types_allowed = True

    @property
    def y_delta(self) -> float:
        return self.y_max - self.y_min

    def _plot_title(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        x_pos = (
            gr.start + (0.01 * gr.length)
            if self.title_location == "left"
            else gr.end - (0.01 * gr.length)
        )
        h_align = "left" if self.title_location == "left" else "right"

        ax.text(
            x_pos,
            self.y_delta * self.title_height,
            self.title,
            horizontalalignment=h_align,
            verticalalignment="top",
            fontdict={
                "size": self.title_size,
                "color": self.title_color,
                "fontname": self.title_font,
                "weight": self.title_weight,
            },
        )

    def _format_scale(self, value: float) -> str:
        if value % 1 == 0:
            return str(int(value))
        return f"{value:.{self.scale_precision}f}"

    def _plot_scale(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        y_min = self._format_scale(self.y_min)
        y_max = self._format_scale(self.y_max)

        x_pos = (
            gr.end - (0.01 * gr.length)
            if self.scale_location == "right"
            else gr.start + (0.01 * gr.length)
        )
        h_align = "right" if self.scale_location == "right" else "left"

        ax.text(
            x_pos,
            self.y_delta * self.scale_height,
            f"[ {y_min} - {y_max} ]",
            horizontalalignment=h_align,
            verticalalignment="top",
            fontdict={
                "size": self.scale_size,
                "color": self.scale_color,
                "fontname": self.scale_font,
                "weight": self.scale_weight,
            },
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> "TrackLabeller":
        if self.plot_title and self.title:
            self._plot_title(ax, gr)
        if self.plot_scale:
            self._plot_scale(ax, gr)

        import plotnado.tracks as pg_tracks

        pg_tracks.clean_axis(ax)
        return self
