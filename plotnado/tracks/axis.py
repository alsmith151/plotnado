"""
Genomic axis track for showing coordinate scale.
"""

import matplotlib.axes
import matplotlib.ticker
import numpy as np
from pydantic import BaseModel

from .region import GenomicRegion
from .base import Track


class GenomicAxisAesthetics(BaseModel):
    """
    Aesthetics configuration for genomic axis.

    Attributes:
        color: Color for axis line and labels
        font_size: Font size for tick labels
        num_ticks: Approximate number of ticks to display
        show_chromosome: Whether to show chromosome name
    """

    color: str = "black"
    font_size: int = 10
    num_ticks: int = 5
    show_chromosome: bool = True
    tick_height: float = 0.1


class GenomicAxis(Track):
    """
    Track showing genomic coordinates as an axis.

    Displays tick marks and labels for the genomic region.
    """

    title: str = ""
    aesthetics: GenomicAxisAesthetics = GenomicAxisAesthetics()
    height: float = 0.3

    def fetch_data(self, gr: GenomicRegion) -> None:
        return None

    @staticmethod
    def _format_position(pos: int) -> str:
        """Format a genomic position for display."""
        if pos >= 1_000_000:
            return f"{pos / 1_000_000:.1f}M"
        elif pos >= 1_000:
            return f"{pos / 1_000:.1f}k"
        return str(pos)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the genomic axis."""
        # Calculate nice tick positions
        num_ticks = self.aesthetics.num_ticks
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

        # Draw axis line
        y_line = 0.7
        ax.plot(
            [gr.start, gr.end],
            [y_line, y_line],
            color=self.aesthetics.color,
            linewidth=1,
        )

        # Draw ticks and labels
        tick_height = self.aesthetics.tick_height
        for pos in tick_positions:
            ax.plot(
                [pos, pos],
                [y_line - tick_height, y_line + tick_height],
                color=self.aesthetics.color,
                linewidth=1,
            )
            ax.text(
                pos,
                y_line - tick_height * 2,
                self._format_position(int(pos)),
                ha="center",
                va="top",
                fontsize=self.aesthetics.font_size,
                color=self.aesthetics.color,
            )

        # Draw chromosome label if enabled
        if self.aesthetics.show_chromosome:
            ax.text(
                gr.start,
                0.95,
                gr.chromosome,
                ha="left",
                va="top",
                fontsize=self.aesthetics.font_size,
                color=self.aesthetics.color,
                fontweight="bold",
            )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)

        # Hide all spines and ticks
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
