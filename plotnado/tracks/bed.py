"""
BED track for displaying genomic intervals.
"""

from pathlib import Path
from typing import List, Union

import matplotlib.axes
import matplotlib.patches
import pandas as pd
from pydantic import BaseModel

from .region import GenomicRegion
from .base import Track
from .utils import clean_axis
from .enums import DisplayMode


class BedAesthetics(BaseModel):
    """
    Aesthetics configuration for BED tracks.

    Attributes:
        color: Fill color for intervals
        edge_color: Edge color for intervals
        alpha: Transparency (0-1)
        interval_height: Height of each interval (0-1)
        display: Display mode (collapsed or expanded)
        max_rows: Maximum number of rows for expanded display
    """

    color: str = "steelblue"
    edge_color: str = "black"
    alpha: float = 0.8
    interval_height: float = 0.6
    display: DisplayMode = DisplayMode.COLLAPSED
    max_rows: int = 5
    show_labels: bool = False
    label_field: str = "name"
    font_size: int = 8


class BedTrack(Track):
    """
    Track for displaying BED intervals.

    Supports BED3, BED4, BED6, and BED12 formats.

    Attributes:
        data: Path to BED file or DataFrame with BED-like columns
        aesthetics: Visual styling configuration
    """

    data: Union[Path, pd.DataFrame, str]
    aesthetics: BedAesthetics = BedAesthetics()
    height: float = 1.0

    class Config:
        arbitrary_types_allowed = True

    def _fetch_from_disk(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch intervals from a BED file."""
        from pybedtools import BedTool

        path = str(self.data)
        bt = BedTool(path)

        # Get intervals overlapping region
        region_str = f"{gr.chromosome}:{gr.start}-{gr.end}"
        try:
            bt_tabix = bt.tabix(force=True)
            intervals = bt_tabix.tabix_intervals(region_str)
            df = intervals.to_dataframe()
        except (OSError, Exception):
            # Fallback: filter in memory
            df = bt.to_dataframe()
            if not df.empty:
                df = df[
                    (df["chrom"] == gr.chromosome)
                    & (df["end"] > gr.start)
                    & (df["start"] < gr.end)
                ]

        return df

    def _fetch_from_df(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch intervals from a DataFrame."""
        df = self.data
        # Handle different column names
        chrom_col = "chrom" if "chrom" in df.columns else "Chromosome"
        start_col = "start" if "start" in df.columns else "Start"
        end_col = "end" if "end" in df.columns else "End"

        mask = (
            (df[chrom_col] == gr.chromosome)
            & (df[end_col] > gr.start)
            & (df[start_col] < gr.end)
        )
        result = df.loc[mask].copy()

        # Normalize column names
        if "Chromosome" in result.columns:
            result = result.rename(
                columns={"Chromosome": "chrom", "Start": "start", "End": "end"}
            )
        return result

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch BED data for the given region."""
        if isinstance(self.data, pd.DataFrame):
            return self._fetch_from_df(gr)
        return self._fetch_from_disk(gr)

    def _allocate_row_index(
        self, row_last_positions: List[int], start_bp: int, end_bp: int
    ) -> int:
        """Allocate a row where the new interval doesn't overlap existing ones."""
        for idx, last_end in enumerate(row_last_positions):
            if last_end < start_bp:
                row_last_positions[idx] = end_bp
                return idx
        row_last_positions.append(end_bp)
        return len(row_last_positions) - 1

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot BED records."""
        data = self.fetch_data(gr)

        if data.empty:
            ax.set_xlim(gr.start, gr.end)
            ax.set_ylim(0, 1)
            clean_axis(ax)
            return

        row_scale = 1.0 / max(1, self.aesthetics.max_rows)
        row_last_positions: List[int] = []

        for row in data.itertuples():
            start = row.start
            end = row.end

            if self.aesthetics.display == DisplayMode.COLLAPSED:
                row_index = 0
            else:
                row_index = self._allocate_row_index(row_last_positions, start, end)

            if row_index >= self.aesthetics.max_rows:
                continue

            ypos = (
                0.5
                if self.aesthetics.display == DisplayMode.COLLAPSED
                else ((row_index + 0.5) * row_scale)
            )

            # Draw interval
            rect = matplotlib.patches.Rectangle(
                (start, ypos - self.aesthetics.interval_height / 2),
                end - start,
                self.aesthetics.interval_height,
                linewidth=1,
                edgecolor=self.aesthetics.edge_color,
                facecolor=self.aesthetics.color,
                alpha=self.aesthetics.alpha,
            )
            ax.add_patch(rect)

            # Draw label if enabled
            if self.aesthetics.show_labels and hasattr(
                row, self.aesthetics.label_field
            ):
                label = getattr(row, self.aesthetics.label_field)
                ax.text(
                    (start + end) / 2,
                    ypos,
                    str(label),
                    ha="center",
                    va="center",
                    fontsize=self.aesthetics.font_size,
                )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
        clean_axis(ax)
