"""
Scale bar track for showing genomic distances.
"""

import matplotlib.axes
import matplotlib.patches
from pydantic import Field

from .region import GenomicRegion
from .base import Track
from .utils import clean_axis, format_distance
from .enums import PlotStyle, Position, TrackType
from .aesthetics import BaseAesthetics
from .registry import registry


class ScaleBarAesthetics(BaseAesthetics):
    """Visual configuration for scale bars.

    Inherits color, alpha, and linewidth from BaseAesthetics.
    """

    style: PlotStyle = Field(default=PlotStyle.STD, description="Scale bar style variant.")
    position: Position = Field(default=Position.LEFT, description="Horizontal anchor of the scale bar.")
    scale_distance: float | None = Field(
        default=None,
        description="Explicit scale bar length in base pairs; auto-selected if omitted.",
    )
    font_size: int = Field(default=8, description="Font size for the scale label.")
    title: str = Field(default="Scale", description="Human-readable name for this style preset.")
    bar_linewidth: float = Field(default=1.8, description="Line width of the main horizontal scale bar.")
    tick_linewidth: float = Field(default=1.4, description="Line width of terminal scale ticks.")
    tick_height: float = Field(default=0.1, description="Half-height of terminal ticks in axis units.")
    label_offset: float = Field(
        default=0.25,
        description="Vertical distance between bar baseline and text label.",
    )


@registry.register(TrackType.SCALEBAR, aliases=["scale"])
class ScaleBar(Track):
    """
    Track showing a scale bar with distance annotation.
    """

    title: str = Field(default="ScaleBar", description="Optional title for the scale bar track.")
    aesthetics: ScaleBarAesthetics = Field(
        default_factory=ScaleBarAesthetics,
        description="Visual style options for the scale bar and label.",
    )
    height: float = Field(default=0.3, description="Relative panel height for this compact track.")

    def fetch_data(self, gr: GenomicRegion) -> None:
        return None

    @staticmethod
    def _get_appropriate_scale(length: int) -> int:
        """Determine an appropriate scale for a genomic region of given length."""
        if length <= 0:
            raise ValueError("Length must be positive")

        target = max(1, int(length * 0.2))
        nice = ScaleBar._round_to_nice_scale(target)
        return min(nice, length)

    @staticmethod
    def _round_to_nice_scale(value: int) -> int:
        if value <= 0:
            return 1

        candidates = []
        magnitude = 1
        while magnitude <= value * 10:
            for base in (1, 2, 5):
                candidates.append(base * magnitude)
            magnitude *= 10

        below_or_equal = [candidate for candidate in candidates if candidate <= value]
        if below_or_equal:
            return max(below_or_equal)
        return min(candidates)

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
            solid_capstyle="butt",
            clip_on=False,
            zorder=5,
        )

        # Plot ticks (as expected by tests)
        ax.plot(
            [x0, x0],
            [y_midpoint - self.tick_height, y_midpoint + self.tick_height],
            color=color,
            linewidth=self.tick_linewidth,
            solid_capstyle="butt",
            clip_on=False,
            zorder=5,
        )
        ax.plot(
            [x1, x1],
            [y_midpoint - self.tick_height, y_midpoint + self.tick_height],
            color=color,
            linewidth=self.tick_linewidth,
            solid_capstyle="butt",
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
