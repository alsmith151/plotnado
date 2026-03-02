"""
Links track for visualizing genomic interactions (arcs).
"""

from pathlib import Path

import matplotlib.axes
import matplotlib.cm
import matplotlib.colors
import matplotlib.patches
import matplotlib.path
import pandas as pd
from pydantic import BaseModel, ConfigDict

from .region import GenomicRegion
from .base import Track
from .utils import clean_axis


class LinksAesthetics(BaseModel):
    """
    Aesthetics configuration for Links tracks.

    Attributes:
        color: Default color for arcs
        edge_color: Edge color for arcs (usually same as color)
        alpha: Transparency (0-1)
        linewidth: Width of arc lines
        cmap: Colormap to use for score-based coloring
        max_height: Maximum height of arcs relative to track height (0-1)
        color_by_score: Whether to use the score column for coloring
    """

    color: str = "steelblue"
    edge_color: str | None = None
    alpha: float = 0.6
    linewidth: float = 1.0
    cmap: str = "viridis"
    max_height: float = 0.8
    color_by_score: bool = False
    min_score: float | None = None
    max_score: float | None = None
    y_baseline: float = 0.1


class LinksTrack(Track):
    """
    Track for displaying genomic links/arcs (e.g., Hi-C loops, ChIA-PET).

    Supports BEDPE-like data: chrom1, start1, end1, chrom2, start2, end2, [score].

    Attributes:
        data: Path to BEDPE file or DataFrame
        aesthetics: Visual styling configuration
    """

    data: Path | pd.DataFrame | str
    aesthetics: LinksAesthetics = LinksAesthetics()
    height: float = 2.0

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch links overlapping the region."""
        if isinstance(self.data, pd.DataFrame):
            df = self.data
        else:
            # For simplicity, assume path points to a file that can be read by pandas
            # In a real implementation, we'd use pybedtools or similar for large files
            df = pd.read_csv(self.data, sep="\t", header=None)

            # Map standard BEDPE columns if no header
            if df.columns.dtype == int:
                cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
                if len(df.columns) > 6:
                    cols.append("score")
                df.columns = cols[: len(df.columns)]

        # Filter for intra-chromosomal links in current region
        mask = (
            (df["chrom1"] == gr.chromosome)
            & (df["chrom2"] == gr.chromosome)
            & (
                ((df["start1"] >= gr.start) & (df["start1"] <= gr.end))
                | ((df["end1"] >= gr.start) & (df["end1"] <= gr.end))
                | ((df["start2"] >= gr.start) & (df["start2"] <= gr.end))
                | ((df["end2"] >= gr.start) & (df["end2"] <= gr.end))
            )
        )
        return df.loc[mask].copy()

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot links as arcs."""
        data = self.fetch_data(gr)

        if data.empty:
            ax.set_xlim(gr.start, gr.end)
            ax.set_ylim(0, 1)
            clean_axis(ax)
            return

        # Setup coloring
        cmap = None
        norm = None
        if self.color_by_score and "score" in data.columns:
            cmap = matplotlib.cm.get_cmap(self.cmap)
            vmin = (
                self.min_score
                if self.min_score is not None
                else data["score"].min()
            )
            vmax = (
                self.max_score
                if self.max_score is not None
                else data["score"].max()
            )
            if vmin == vmax:
                vmax += 1
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

        # Baseline for arcs
        y_baseline = self.y_baseline

        # Calculate max distance for height scaling
        # Use region width as a proxy for normalization if needed, or constant factor
        region_width = gr.end - gr.start

        for row in data.itertuples():
            # Get centers of both anchors
            x1 = (row.start1 + row.end1) / 2
            x2 = (row.start2 + row.end2) / 2

            # Distance
            dist = abs(x2 - x1)
            if dist == 0:
                continue

            # Height proportional to distance, capped at max_height
            # Scaling factor can be adjusted. Let's try quadratic or sqrt for better visibility of short/long links
            h = min(self.max_height, (dist / region_width) * 2)

            # Calculate Bezier control point
            xm = (x1 + x2) / 2
            ym = y_baseline + h

            # Color
            color = self.color
            if cmap and norm:
                color = cmap(norm(row.score))

            # Draw quadratic Bezier arc
            path_data = [
                (matplotlib.path.Path.MOVETO, (x1, y_baseline)),
                (matplotlib.path.Path.CURVE3, (xm, ym)),
                (matplotlib.path.Path.CURVE3, (x2, y_baseline)),
            ]

            codes, verts = zip(*path_data)
            path = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(
                path,
                facecolor="none",
                edgecolor=color,
                alpha=self.alpha,
                linewidth=self.linewidth,
                zorder=1,
            )
            ax.add_patch(patch)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)  # Internal coordinate system for track height
        clean_axis(ax)
