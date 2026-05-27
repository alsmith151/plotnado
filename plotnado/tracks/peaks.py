"""
NarrowPeak track for visualizing ChIP-seq/ATAC-seq peaks.
"""

import matplotlib.axes
import matplotlib.cm
import matplotlib.colors
import matplotlib.patches
from pydantic import Field

from .bed import BedTrack, BedAesthetics
from .base import TrackLabeller
from .region import GenomicRegion
from .utils import clean_axis
from .enums import DisplayMode, NarrowPeakColorBy, TrackType
from .registry import registry


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

    color: str = Field(default="#d95f02", description="Fallback color for peaks.")
    interval_height: float = Field(
        default=0.28,
        description="Rectangle height in normalized track coordinates for peak intervals.",
    )
    color_by: NarrowPeakColorBy | None = Field(
        default=NarrowPeakColorBy.SIGNAL_VALUE,
        description="Optional narrowPeak field used for colormap-based coloring.",
    )
    cmap: str = Field(default="Oranges", description="Colormap used when color_by is enabled.")
    min_score: float | None = Field(
        default=None,
        description="Optional lower score bound for colormap normalization.",
    )
    max_score: float | None = Field(
        default=None,
        description="Optional upper score bound for colormap normalization.",
    )
    show_summit: bool = Field(default=True, description="Draw summit marker when peak offset is available.")
    summit_color: str = Field(default="black", description="Color of summit marker lines.")
    summit_width: float = Field(default=0.8, description="Line width of summit marker lines.")


@registry.register(TrackType.NARROWPEAK)
class NarrowPeakTrack(BedTrack):
    """
    Track for displaying NarrowPeak intervals.

    Attributes:
        data: Path to narrowPeak file or DataFrame
        aesthetics: Visual styling configuration
    """

    aesthetics: NarrowPeakAesthetics = Field(
        default_factory=NarrowPeakAesthetics,
        description="NarrowPeak-specific visual styling options.",
    )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot narrowPeak records."""
        data = self.fetch_data(gr)
        peak_label = self.label.model_copy(update={"plot_scale": False})

        if data.empty:
            ax.set_xlim(gr.start, gr.end)
            ax.set_ylim(0, 1)
            if self.label.plot_title or self.label.plot_scale:
                TrackLabeller.from_config(
                    peak_label,
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

        # Setup colormap if needed
        cmap = None
        norm = None
        if self.color_by and self.color_by in data.columns:
            values = data[self.color_by]
            cmap = matplotlib.cm.get_cmap(self.cmap)
            vmin = (
                self.min_score
                if self.min_score is not None
                else values.min()
            )
            vmax = (
                self.max_score
                if self.max_score is not None
                else values.max()
            )
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

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

            # Determine color
            current_color = self.color
            if cmap and norm and hasattr(row, self.color_by):
                val = getattr(row, self.color_by)
                current_color = cmap(norm(val))

            # Draw interval
            rect = matplotlib.patches.Rectangle(
                (start, ypos - self.interval_height / 2),
                end - start,
                self.interval_height,
                linewidth=self.rect_linewidth if self.draw_edges else 0,
                edgecolor=self.edge_color if self.draw_edges else "none",
                facecolor=current_color,
                alpha=self.alpha,
                zorder=1,
            )
            ax.add_patch(rect)

            # Draw summit if enabled
            if self.show_summit and hasattr(row, "peak") and row.peak != -1:
                # peak is 0-based offset from start
                summit_pos = start + row.peak
                ax.plot(
                    [summit_pos, summit_pos],
                    [
                        ypos - self.interval_height / 2,
                        ypos + self.interval_height / 2,
                    ],
                    color=self.summit_color,
                    linewidth=self.summit_width,
                    zorder=2,
                )

            # Draw label if enabled
            if self.show_labels and hasattr(row, self.label_field):
                label = getattr(row, self.label_field)
                ax.text(
                    (start + end) / 2,
                    ypos,
                    str(label),
                    ha="center",
                    va="center",
                    fontsize=self.font_size,
                    zorder=3,
                )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
        if self.label.plot_title or self.label.plot_scale:
            TrackLabeller.from_config(
                peak_label,
                gr,
                0,
                1,
                title=self.title or "",
                title_color=self.color,
            ).plot(ax, gr)
        else:
            clean_axis(ax)
