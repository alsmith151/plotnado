"""
BigWig track for signal visualization.
"""

from pathlib import Path
from typing import Literal, Optional, Union

import matplotlib.axes
import matplotlib.patches
import numpy as np
import pandas as pd
import pybigtools
from pydantic import BaseModel

from .region import GenomicRegion
from .base import Track, TrackLabeller
from .utils import clean_axis
from .schemas import BedgraphDataFrame


class BigwigAesthetics(BaseModel):
    """
    Aesthetics configuration for BigWig tracks.

    Attributes:
        style: Plot style ('std', 'line', 'scatter', 'heatmap')
    """

    style: Literal["std", "scatter", "line", "heatmap"] = "std"
    color: str = "black"
    fill: bool = True
    alpha: float = 1.0
    linewidth: float = 1.0
    scatter_point_size: float = 1.0

    min_value: Optional[float] = None
    max_value: Optional[float] = None

    plot_title: bool = True
    title_location: Literal["left", "right"] = "left"
    title_height: float = 0.5  # Add this field expected by some old tests

    plot_scale: bool = True
    scale_location: Literal["left", "right"] = "left"
    scale_height: float = 0.5
    title: str = "BigWig"


class BigWigTrack(Track):
    """
    Track for displaying BigWig signal data.

    Attributes:
        data: Path to BigWig file or DataFrame with bedgraph-like data
        aesthetics: Visual styling configuration
    """

    data: Optional[Union[Path, pd.DataFrame, str]] = None
    aesthetics: BigwigAesthetics = BigwigAesthetics()

    y_min: Optional[float] = None
    y_max: Optional[float] = None

    def _plot_stairs(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, values: pd.DataFrame
    ) -> None:
        """Plot as stairs."""
        edges = np.linspace(gr.start, gr.end, values.shape[0] + 1)
        ax.stairs(
            values=values["value"],
            edges=edges,
            linewidth=self.aesthetics.linewidth,
            color=self.aesthetics.color,
            alpha=self.aesthetics.alpha,
            fill=self.aesthetics.fill,
        )

    def _fetch_from_disk(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch data from a BigWig file."""
        path = str(self.data)
        bw = pybigtools.open(path)
        records = list(bw.records(gr.chromosome, gr.start, gr.end))

        if not records:
            return pd.DataFrame(columns=["chrom", "start", "end", "value"])

        df = pd.DataFrame(records, columns=["start", "end", "value"])
        df["chrom"] = gr.chromosome
        return BedgraphDataFrame(df[["start", "end", "value", "chrom"]])

    def _fetch_from_df(self, gr: GenomicRegion) -> BedgraphDataFrame:
        """Fetch data from a DataFrame."""
        df = self.data
        mask = (
            (df["chrom"] == gr.chromosome)
            & (df["end"] > gr.start)
            & (df["start"] < gr.end)
        )
        return BedgraphDataFrame(
            df.loc[mask, ["start", "end", "value", "chrom"]].copy()
        )

    def fetch_data(self, gr: GenomicRegion) -> BedgraphDataFrame:
        """Fetch data for the given genomic region."""
        if isinstance(self.data, pd.DataFrame):
            df = self._fetch_from_df(gr)
        else:
            df = self._fetch_from_disk(gr)
        return BedgraphDataFrame(df)

    def _plot_fill(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, data: pd.DataFrame
    ) -> None:
        """Plot as filled area."""
        if data.empty:
            return

        import numpy as np

        # Create step plot using start/end positions
        x = []
        y = []
        for _, row in data.iterrows():
            x.extend([row["start"], row["end"]])
            y.extend([row["value"], row["value"]])

        ax.fill_between(
            x,
            y,
            alpha=self.aesthetics.alpha,
            color=self.aesthetics.color,
            linewidth=self.aesthetics.linewidth,
            step="pre" if self.aesthetics.style == "std" else None,
        )

    def _plot_line(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, data: pd.DataFrame
    ) -> None:
        """Plot as line."""
        if data.empty:
            return

        x = (data["start"] + data["end"]) / 2
        y = data["value"]
        ax.plot(
            x,
            y,
            color=self.aesthetics.color,
            alpha=self.aesthetics.alpha,
            linewidth=self.aesthetics.linewidth,
        )

    def _plot_scatter(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, data: pd.DataFrame
    ) -> None:
        """Plot as scatter."""
        if data.empty:
            return

        x = (data["start"] + data["end"]) / 2
        y = data["value"]
        ax.scatter(
            x,
            y,
            color=self.aesthetics.color,
            alpha=self.aesthetics.alpha,
            s=self.aesthetics.scatter_point_size,
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the BigWig track."""
        data = self.fetch_data(gr)

        # Plot based on style
        if self.aesthetics.style == "scatter":
            self._plot_scatter(ax, gr, data)
        elif self.aesthetics.style == "line":
            self._plot_line(ax, gr, data)
        elif self.aesthetics.style == "heatmap":
            pass  # Heatmap not implemented yet

        if self.aesthetics.style == "std":
            self._plot_stairs(ax, gr, data)

        # Set axis limits
        ax.set_xlim(gr.start, gr.end)

        if data.empty:
            self.y_min = 0
            self.y_max = 1
        else:
            self.y_min = (
                self.aesthetics.min_value
                if self.aesthetics.min_value is not None
                else data["value"].min()
            )
            self.y_max = (
                self.aesthetics.max_value
                if self.aesthetics.max_value is not None
                else data["value"].max()
            )

        # Ensure some range
        if self.y_min == self.y_max:
            self.y_max = self.y_min + 1

        ax.set_ylim(ymin=self.y_min, ymax=self.y_max)

        # Add labels
        if self.aesthetics.plot_title or self.aesthetics.plot_scale:
            labeller = TrackLabeller(
                gr=gr,
                y_min=self.y_min,
                y_max=self.y_max,
                plot_title=self.aesthetics.plot_title,
                plot_scale=self.aesthetics.plot_scale,
                title=self.title or "",
                scale_min=self.y_min,
                scale_max=self.y_max,
                title_location=self.aesthetics.title_location,
                scale_location=self.aesthetics.scale_location,
            )
            labeller.plot(ax, gr)
        else:
            clean_axis(ax)
