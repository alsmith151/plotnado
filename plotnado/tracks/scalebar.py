"""
Scale bar track for showing genomic distances.
"""

from typing import Literal

import matplotlib.axes
import matplotlib.patches
from pydantic import BaseModel

from .region import GenomicRegion
from .base import Track
from .utils import clean_axis, format_distance
from .enums import Position


class ScaleBarAesthetics(BaseModel):
    """Visual configuration for scale bars."""

    style: Literal["std"] = "std"
    color: str = "#333333"  # Dark gray
    position: Position | Literal["left", "right", "center"] = Position.LEFT
    scale_distance: float | None = None
    font_size: int = 8
    title: str = "Scale"
    bar_linewidth: float = 3.0
    tick_linewidth: float = 2.0
    tick_height: float = 0.1
    label_offset: float = 0.25


class ScaleBar(Track):
    """
    Track showing a scale bar with distance annotation.
    """

    title: str = "ScaleBar"
    aesthetics: ScaleBarAesthetics = ScaleBarAesthetics()
    height: float = 0.3

    def fetch_data(self, gr: GenomicRegion) -> None:
        return None

    @staticmethod
    def _get_appropriate_scale(length: int) -> int:
        """Determine an appropriate scale for a genomic region of given length."""
        if length <= 0:
            raise ValueError("Length must be positive")

        if length <= 100:
            return 10
        elif length <= 1000:
            return 200
        elif length <= 10000:
            return 2000
        elif length <= 100000:
            return 20000

        # Default to 20% of length for larger regions
        return int(length * 0.2)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the scale bar."""
        position = self.position
        y_midpoint = 0.5
        color = self.color

        scale_distance = self.scale_distance or self._get_appropriate_scale(gr.length)

        # Determine x start and end based on position
        if position == "left":
            x0 = gr.start
            x1 = x0 + scale_distance
        elif position == "right":
            x0 = gr.end - scale_distance
            x1 = gr.end
        elif position == "center":
            x0 = gr.center - (scale_distance / 2)
            x1 = gr.center + (scale_distance / 2)
        else:
            raise ValueError('Position can only be "left", "right" or "center"')

        # Plot scale bar (thicker for visibility)
        ax.plot(
            [x0, x1],
            [y_midpoint, y_midpoint],
            color=color,
            linewidth=self.bar_linewidth,
            clip_on=False,
            zorder=5,
        )

        # Plot ticks (as expected by tests)
        ax.plot(
            [x0, x0],
            [y_midpoint - self.tick_height, y_midpoint + self.tick_height],
            color=color,
            linewidth=self.tick_linewidth,
            clip_on=False,
            zorder=5,
        )
        ax.plot(
            [x1, x1],
            [y_midpoint - self.tick_height, y_midpoint + self.tick_height],
            color=color,
            linewidth=self.tick_linewidth,
            clip_on=False,
            zorder=5,
        )

        # Add label
        scale_label = format_distance(scale_distance)
        ax.text(
            (x0 + x1) / 2,
            y_midpoint - self.label_offset,
            scale_label,
            ha="center",
            va="center",
            fontsize=self.font_size,
        )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
        clean_axis(ax)
