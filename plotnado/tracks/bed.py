"""
BED track for displaying genomic intervals.
"""

from pathlib import Path
from typing import Any

import matplotlib.axes
import matplotlib.patches
import pandas as pd
from pydantic import ConfigDict, Field, BaseModel

from .region import GenomicRegion
from .base import Track, TrackLabeller
from .utils import clean_axis, read_bed_regions
from .enums import DisplayMode, TrackType
from .aesthetics import BaseAesthetics
from .registry import registry


class BedAesthetics(BaseAesthetics):
    """
    Aesthetics configuration for BED tracks.

    Inherits color, alpha, and linewidth from BaseAesthetics.

    Attributes:
        edge_color: Edge color for intervals
        interval_height: Height of each interval (0-1)
        display: Display mode (collapsed or expanded)
        max_rows: Maximum number of rows for expanded display
    """

    edge_color: str = Field(default="black", description="Stroke color for interval borders.")
    interval_height: float = Field(
        default=0.45,
        description="Rectangle height in normalized track coordinates.",
    )
    display: DisplayMode = Field(
        default=DisplayMode.COLLAPSED,
        description="Collapsed draws all intervals on one row; expanded stacks overlaps.",
    )
    max_rows: int = Field(default=5, description="Maximum stacked rows when display is expanded.")
    show_labels: bool = Field(default=False, description="Show text labels for intervals.")
    label_field: str = Field(
        default="name",
        description="Column name used to populate interval labels.",
    )
    font_size: int = Field(default=8, description="Font size for interval labels.")
    rect_linewidth: float = Field(default=0.7, description="Border line width for interval rectangles.")
    draw_edges: bool = Field(default=True, description="Draw rectangle borders for intervals.")


@registry.register(TrackType.BED, aliases=["annotation", "unknown"])
class BedTrack(Track):
    """
    Track for displaying BED intervals.

    Supports BED3, BED4, BED6, and BED12 formats.

    Attributes:
        data: Path to BED file or DataFrame with BED-like columns
        aesthetics: Visual styling configuration
    """

    data: Path | pd.DataFrame | str | Any | None = Field(
        default=None,
        description="BED/BigBed path or in-memory interval table-like object.",
    )
    aesthetics: BedAesthetics = Field(
        default_factory=BedAesthetics,
        description="Visual styling options for interval rendering.",
    )
    height: float = Field(default=1.0, description="Relative panel height for this track.")

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def _fetch_from_disk(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch intervals from a BED/BigBed file."""
        return read_bed_regions(str(self.data), gr.chromosome, gr.start, gr.end)

    def _fetch_from_memory(self, gr: GenomicRegion, data: Any) -> pd.DataFrame:
        """Fetch intervals from in-memory DataFrame/PyRanges-like object."""
        if hasattr(data, "df"):
            df = data.df
        elif hasattr(data, "as_df"):
            df = data.as_df()
        else:
            df = data

        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)

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
        if self.data is None:
            return pd.DataFrame(columns=["chrom", "start", "end"])

        if isinstance(self.data, pd.DataFrame) or hasattr(self.data, "df"):
            return self._fetch_from_memory(gr, self.data)
        return self._fetch_from_disk(gr)

    def _allocate_row_index(
        self, row_last_positions: list[int], start_bp: int, end_bp: int
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
        bed_label = self.label.model_copy(update={"plot_scale": False})

        if data.empty:
            ax.set_xlim(gr.start, gr.end)
            ax.set_ylim(0, 1)
            if self.label.plot_title or self.label.plot_scale:
                TrackLabeller.from_config(
                    bed_label,
                    gr,
                    0,
                    1,
                    title=self.title or "",
                    title_color=self.color,
                ).plot(ax, gr)
            else:
                clean_axis(ax)
            return

        row_scale = 1.0 / max(1, self.max_rows)
        row_last_positions: list[int] = []

        for row in data.itertuples():
            start = row.start
            end = row.end

            if self.display == DisplayMode.COLLAPSED:
                row_index = 0
            else:
                row_index = self._allocate_row_index(row_last_positions, start, end)

            if row_index >= self.max_rows:
                continue

            ypos = (
                0.5
                if self.display == DisplayMode.COLLAPSED
                else ((row_index + 0.5) * row_scale)
            )

            # Draw interval
            rect = matplotlib.patches.Rectangle(
                (start, ypos - self.interval_height / 2),
                end - start,
                self.interval_height,
                linewidth=self.rect_linewidth if self.draw_edges else 0,
                edgecolor=self.edge_color if self.draw_edges else "none",
                facecolor=self.color,
                alpha=self.alpha,
            )
            ax.add_patch(rect)

            # Draw label if enabled
            if self.show_labels and hasattr(row, self.label_field):
                label = getattr(row, self.label_field)
                # Position label above the peak, within track bounds
                label_ypos = ypos + self.interval_height / 2 + 0.05
                ax.text(
                    (start + end) / 2,
                    label_ypos,
                    str(label),
                    ha="center",
                    va="bottom",
                    fontsize=self.font_size,
                    clip_on=True,  # Clip text that extends outside axis
                )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
        if self.label.plot_title or self.label.plot_scale:
            TrackLabeller.from_config(
                bed_label,
                gr,
                0,
                1,
                title=self.title or "",
                title_color=self.color,
            ).plot(ax, gr)
        else:
            clean_axis(ax)
