"""
Spacer track for adding vertical space between tracks.
"""

import matplotlib.axes
from .base import GenomicRegion, Track
from .utils import clean_axis


class Spacer(Track):
    """
    Empty spacer track for adding vertical space.
    """

    title: str = ""
    height: float = 0.5

    def fetch_data(self, gr: GenomicRegion) -> None:
        return None

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot empty spacer."""
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
        clean_axis(ax)
