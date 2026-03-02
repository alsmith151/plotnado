"""
Highlights track for marking genomic regions.
"""

from pathlib import Path

import matplotlib.axes
import matplotlib.patches
import pandas as pd
from pydantic import BaseModel, ConfigDict

from .region import GenomicRegion
from .base import Track
from .utils import read_bed_regions


class HighlightsAesthetics(BaseModel):
    """
    Aesthetics configuration for highlight regions.

    Attributes:
        color: Fill color for highlighted regions
        alpha: Transparency (0-1)
        edge_color: Edge color (None for no edge)
    """

    color: str = "yellow"
    alpha: float = 0.3
    edge_color: str | None = None
    linewidth: float = 1.0


class HighlightsFromFile(Track):
    """
    Track for highlighting genomic regions from a BED file.

    This creates transparent overlays on the specified regions.
    Useful for marking regions of interest, peaks, or annotations.

    Attributes:
        data: Path to BED file or DataFrame with regions to highlight
        aesthetics: Visual styling configuration
    """

    data: Path | pd.DataFrame | str
    aesthetics: HighlightsAesthetics = HighlightsAesthetics()
    height: float = 0.0  # Zero height means it spans all tracks

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def _fetch_from_disk(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch highlight regions from BED/BigBed."""
        return read_bed_regions(str(self.data), gr.chromosome, gr.start, gr.end)

    def _fetch_from_df(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch highlight regions from a DataFrame."""
        df = self.data
        chrom_col = "chrom" if "chrom" in df.columns else "Chromosome"
        start_col = "start" if "start" in df.columns else "Start"
        end_col = "end" if "end" in df.columns else "End"

        mask = (
            (df[chrom_col] == gr.chromosome)
            & (df[end_col] > gr.start)
            & (df[start_col] < gr.end)
        )
        result = df.loc[mask].copy()

        if "Chromosome" in result.columns:
            result = result.rename(
                columns={"Chromosome": "chrom", "Start": "start", "End": "end"}
            )
        return result

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch regions to highlight."""
        if isinstance(self.data, pd.DataFrame):
            return self._fetch_from_df(gr)
        return self._fetch_from_disk(gr)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot genomic highlights."""
        data = self.fetch_data(gr)

        if data.empty:
            return

        y_min, y_max = ax.get_ylim()
        if y_min == 0 and y_max == 1:
            # Default range if not set
            y_min, y_max = 0, 1

        for row in data.itertuples():
            start = max(row.start, gr.start)
            end = min(row.end, gr.end)

            rect = matplotlib.patches.Rectangle(
                (start, y_min),
                end - start,
                y_max - y_min,
                facecolor=self.color,
                alpha=self.alpha,
                edgecolor=self.edge_color,
                linewidth=0 if self.edge_color is None else self.linewidth,
            )
            ax.add_patch(rect)

    def plot_on_axes(self, gr: GenomicRegion, axes: list) -> None:
        """Plot highlights on multiple axes (for use in Figure)."""
        for ax in axes:
            self.plot(ax, gr)
