"""
Annotation tracks for horizontal and vertical reference lines.
"""

from typing import Any

import matplotlib.axes
from pydantic import BaseModel

from .region import GenomicRegion
from .base import Track


class AnnotationAesthetics(BaseModel):
    """
    Aesthetics configuration for annotation lines.

    Attributes:
        color: Line color
        linestyle: Line style ('-', '--', ':', etc.)
        linewidth: Line width
        alpha: Transparency (0-1)
        zorder: Plotting order
    """

    color: str = "red"
    linestyle: str = "--"
    linewidth: float = 1.0
    alpha: float = 0.8
    zorder: int = 10


class HLineTrack(Track):
    """
    Track for drawing a horizontal reference line.

    Attributes:
        y_value: Y-coordinate for the line
        aesthetics: Visual configuration
    """

    y_value: float
    aesthetics: AnnotationAesthetics = AnnotationAesthetics()
    height: float = 0.0  # Usually overlaid, so zero standalone height

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot horizontal line."""
        ax.axhline(
            y=self.y_value,
            color=self.color,
            linestyle=self.linestyle,
            linewidth=self.linewidth,
            alpha=self.alpha,
            zorder=self.zorder,
        )


class VLineTrack(Track):
    """
    Track for drawing a vertical reference line at a genomic position.

    Attributes:
        x_position: Genomic coordinate for the line
        aesthetics: Visual configuration
    """

    x_position: int | str
    aesthetics: AnnotationAesthetics = AnnotationAesthetics()
    height: float = 0.0

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot vertical line."""
        from .utils import parse_genomic_value

        x = parse_genomic_value(self.x_position)
        if gr.start <= x <= gr.end:
            ax.axvline(
                x=x,
                color=self.color,
                linestyle=self.linestyle,
                linewidth=self.linewidth,
                alpha=self.alpha,
                zorder=self.zorder,
            )

    def plot_on_axes(self, gr: GenomicRegion, axes: list) -> None:
        """Plot vertical line on multiple axes (global marker)."""
        for ax in axes:
            self.plot(ax, gr)
