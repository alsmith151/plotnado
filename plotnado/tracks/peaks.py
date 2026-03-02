"""
NarrowPeak track for visualizing ChIP-seq/ATAC-seq peaks.
"""

from typing import List, Optional, Literal

import matplotlib.axes
import matplotlib.cm
import matplotlib.colors
import matplotlib.patches
import pandas as pd

from .bed import BedTrack, BedAesthetics
from .region import GenomicRegion
from .utils import clean_axis
from .enums import DisplayMode


class NarrowPeakAesthetics(BedAesthetics):
    """
    Aesthetics configuration for NarrowPeak tracks.

    Attributes:
        color_by: Field to use for coloring ('score', 'signalValue', 'pValue', 'qValue', None)
        cmap: Colormap to use when color_by is set
        min_score: Minimum score for color scaling (default: auto)
        max_score: Maximum score for color scaling (default: auto)
        show_summit: Whether to draw a vertical line at the peak summit
    """

    color: str = "#d95f02"  # Default orange for peaks
    color_by: Optional[Literal["score", "signalValue", "pValue", "qValue"]] = (
        "signalValue"
    )
    cmap: str = "Oranges"
    min_score: Optional[float] = None
    max_score: Optional[float] = None
    show_summit: bool = True
    summit_color: str = "black"
    summit_width: float = 1.0


class NarrowPeakTrack(BedTrack):
    """
    Track for displaying NarrowPeak intervals.

    Attributes:
        data: Path to narrowPeak file or DataFrame
        aesthetics: Visual styling configuration
    """

    aesthetics: NarrowPeakAesthetics = NarrowPeakAesthetics()

    def _fetch_from_disk(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch intervals from a narrowPeak file."""
        # Reuse BedTrack logic but ensure columns are named correctly for narrowPeak
        df = super()._fetch_from_disk(gr)

        # If columns are unnamed/standard BED names, try to map them to narrowPeak
        # Standard pybedtools to_dataframe() headers for narrowPeak might vary
        # Usually: chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak

        expected_cols = [
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "signalValue",
            "pValue",
            "qValue",
            "peak",
        ]

        if df.shape[1] == 10:
            # Assuming standard narrowPeak format if exact column count matches
            # but pybedtools usually names the first 6 correctly
            current_cols = list(df.columns)
            mapping = {}
            for i, col in enumerate(expected_cols):
                if i < len(current_cols):
                    mapping[current_cols[i]] = col

            # Only rename if numeric indices or mismatch
            # But proceed with caution.
            pass

        return df

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot narrowPeak records."""
        data = self.fetch_data(gr)

        if data.empty:
            ax.set_xlim(gr.start, gr.end)
            ax.set_ylim(0, 1)
            clean_axis(ax)
            return

        row_scale = 1.0 / max(1, self.aesthetics.max_rows)
        row_last_positions: List[int] = []

        # Setup colormap if needed
        cmap = None
        norm = None
        if self.aesthetics.color_by and self.aesthetics.color_by in data.columns:
            values = data[self.aesthetics.color_by]
            cmap = matplotlib.cm.get_cmap(self.aesthetics.cmap)
            vmin = (
                self.aesthetics.min_score
                if self.aesthetics.min_score is not None
                else values.min()
            )
            vmax = (
                self.aesthetics.max_score
                if self.aesthetics.max_score is not None
                else values.max()
            )
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

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

            # Determine color
            current_color = self.aesthetics.color
            if cmap and norm and hasattr(row, self.aesthetics.color_by):
                val = getattr(row, self.aesthetics.color_by)
                current_color = cmap(norm(val))

            # Draw interval
            rect = matplotlib.patches.Rectangle(
                (start, ypos - self.aesthetics.interval_height / 2),
                end - start,
                self.aesthetics.interval_height,
                linewidth=1,
                edgecolor=self.aesthetics.edge_color,
                facecolor=current_color,
                alpha=self.aesthetics.alpha,
                zorder=1,
            )
            ax.add_patch(rect)

            # Draw summit if enabled
            if self.aesthetics.show_summit and hasattr(row, "peak") and row.peak != -1:
                # peak is 0-based offset from start
                summit_pos = start + row.peak
                ax.plot(
                    [summit_pos, summit_pos],
                    [
                        ypos - self.aesthetics.interval_height / 2,
                        ypos + self.aesthetics.interval_height / 2,
                    ],
                    color=self.aesthetics.summit_color,
                    linewidth=self.aesthetics.summit_width,
                    zorder=2,
                )

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
                    zorder=3,
                )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
        clean_axis(ax)
