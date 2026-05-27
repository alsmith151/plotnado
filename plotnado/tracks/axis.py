"""
Genomic axis track for showing coordinate scale.
"""

import matplotlib.axes
import matplotlib.ticker
import numpy as np
from pydantic import Field

from .region import GenomicRegion
from .base import Track
from .enums import FontWeight, TrackType
from .aesthetics import BaseAesthetics
from .registry import registry


class GenomicAxisAesthetics(BaseAesthetics):
    """
    Aesthetics configuration for genomic axis.

    Inherits color, alpha, and linewidth from BaseAesthetics.

    Attributes:
        font_size: Font size for tick labels
        num_ticks: Approximate number of ticks to display
        show_chromosome: Whether to show chromosome name
    """
    font_size: int = Field(default=9, description="Font size for tick and chromosome labels.")
    num_ticks: int = Field(default=5, description="Target number of tick marks across the region.")
    show_chromosome: bool = Field(default=True, description="Render chromosome name label near the axis.")
    use_human_readable_labels: bool = Field(
        default=False,
        description="Format genomic coordinates using k/M suffixes instead of raw integers.",
    )
    tick_height: float = Field(default=0.15, description="Tick length drawn downward from axis baseline.")
    axis_linewidth: float = Field(default=1.1, description="Line width of the horizontal axis baseline.")
    tick_color: str = Field(default="#333333", description="Color for tick marks and tick labels.")
    tick_linewidth: float = Field(default=0.9, description="Line width of tick marks.")
    chromosome_fontweight: FontWeight = Field(
        default=FontWeight.BOLD,
        description="Font weight of the chromosome text label.",
    )


@registry.register(TrackType.AXIS)
class GenomicAxis(Track):
    """
    Track showing genomic coordinates as an axis.

    Displays tick marks and labels for the genomic region.
    """

    title: str = Field(default="", description="Optional title shown for the axis track.")
    aesthetics: GenomicAxisAesthetics = Field(
        default_factory=GenomicAxisAesthetics,
        description="Visual style options for axis baseline, ticks, and labels.",
    )
    show_chromosome: bool = Field(
        default=False,
        description="Whether to draw the chromosome name on this axis track.",
    )
    height: float = Field(default=0.2, description="Relative panel height for this compact track.")

    def fetch_data(self, gr: GenomicRegion) -> None:
        return None

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the genomic axis."""
        # Calculate nice tick positions
        num_ticks = self.num_ticks
        tick_positions = np.linspace(gr.start, gr.end, num_ticks)

        # Round to nice numbers
        step = (gr.end - gr.start) / (num_ticks - 1)
        magnitude = 10 ** int(np.floor(np.log10(step)))
        nice_step = round(step / magnitude) * magnitude
        first_tick = np.ceil(gr.start / nice_step) * nice_step
        tick_positions = np.arange(first_tick, gr.end, nice_step)

        # Ensure we have at least start and end
        if len(tick_positions) == 0:
            tick_positions = np.array([gr.start, gr.end])

        # Draw axis line at the bottom
        y_line = 0.3
        ax.plot(
            [gr.start, gr.end],
            [y_line, y_line],
            color=self.color,
            linewidth=self.axis_linewidth,
            zorder=2,
        )

        # Draw ticks downward (genome browser convention)
        tick_height = self.tick_height
        for pos in tick_positions:
            tick_label = (
                f"{int(pos):,}"
                if not self.use_human_readable_labels
                else self._format_human_readable_coordinate(int(pos))
            )
            ax.plot(
                [pos, pos],
                [y_line, y_line - tick_height],  # Downward ticks
                color=self.tick_color,
                linewidth=self.tick_linewidth,
                zorder=2,
            )
            ax.text(
                pos,
                y_line - tick_height - 0.05,  # Position below tick
                tick_label,
                ha="center",
                va="top",
                fontsize=self.font_size,
                color=self.tick_color,
                zorder=3,
            )

        # Draw chromosome label if enabled
        if self.show_chromosome:
            ax.text(
                gr.start,
                0.95,
                gr.chromosome,
                ha="left",
                va="top",
                fontsize=self.font_size,
                color=self.color,
                fontweight=self.chromosome_fontweight,
            )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)

        # Hide all spines and ticks
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    @staticmethod
    def _format_human_readable_coordinate(pos: int) -> str:
        if pos >= 1_000_000:
            return f"{(pos / 1_000_000):.2f}".rstrip("0").rstrip(".") + "M"
        if pos >= 1_000:
            return f"{(pos / 1_000):.1f}".rstrip("0").rstrip(".") + "k"
        return str(pos)
