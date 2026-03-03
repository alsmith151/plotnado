"""
BigWig track for signal visualization.
"""

from pathlib import Path

import matplotlib.axes
import matplotlib.patches
import numpy as np
import pandas as pd
import pybigtools
from pydantic import BaseModel
from pydantic import ConfigDict

from .region import GenomicRegion
from .base import Track, TrackLabeller
from .utils import clean_axis
from .schemas import BedgraphDataFrame
from .enums import PlotStyle


class BigwigAesthetics(BaseModel):
    """
    Aesthetics configuration for BigWig tracks.

    Attributes:
        style: Plot style ('std', 'line', 'scatter', 'heatmap')
    """

    style: PlotStyle = PlotStyle.FILL
    color: str = "#2171b5"  # Genome browser blue
    fill: bool = True
    alpha: float = 0.85
    linewidth: float = 1.0
    scatter_point_size: float = 1.0

    min_value: float | None = None
    max_value: float | None = None

    model_config = ConfigDict(use_enum_values=True)


class BigWigTrack(Track):
    """
    Track for displaying BigWig signal data.

    Attributes:
        data: Path to BigWig file or DataFrame with bedgraph-like data
        aesthetics: Visual styling configuration
    """

    data: Path | pd.DataFrame | str | None = None
    aesthetics: BigwigAesthetics = BigwigAesthetics()

    y_min: float | None = None
    y_max: float | None = None

    def _plot_stairs(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, values: pd.DataFrame
    ) -> None:
        """Plot as stairs."""
        edges = np.linspace(gr.start, gr.end, values.shape[0] + 1)
        ax.stairs(
            values=values["value"],
            edges=edges,
            linewidth=self.linewidth,
            color=self.color,
            alpha=self.alpha,
            fill=self.fill,
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
            alpha=self.alpha,
            color=self.color,
            linewidth=0,
            step="post" if self.style in ["std", "fill"] else None,
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
            color=self.color,
            alpha=self.alpha,
            linewidth=self.linewidth,
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
            color=self.color,
            alpha=self.alpha,
            s=self.scatter_point_size,
        )

    def _plot_fragment(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, data: pd.DataFrame
    ) -> None:
        """Plot as explicit fragment blocks using true start/end intervals."""
        if data.empty:
            return

        starts = data["start"].to_numpy(dtype=float)
        widths = (data["end"] - data["start"]).to_numpy(dtype=float)
        values = data["value"].to_numpy(dtype=float)
        ax.bar(
            starts,
            values,
            width=widths,
            align="edge",
            bottom=0,
            color=self.color,
            edgecolor=self.color,
            linewidth=max(self.linewidth, 0.4),
            alpha=self.alpha,
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the BigWig track."""
        data = self.fetch_data(gr)

        # Plot based on style
        if self.style == "scatter":
            self._plot_scatter(ax, gr, data)
        elif self.style == "line":
            self._plot_line(ax, gr, data)
        elif self.style == "fragment":
            self._plot_fragment(ax, gr, data)
        elif self.style in ["std", "fill", "heatmap"]:
            self._plot_fill(ax, gr, data)

        # Set axis limits
        ax.set_xlim(gr.start, gr.end)

        if data.empty:
            self.y_min = 0  # Always start at 0 for genome browser style
            self.y_max = 1
        else:
            # Always use 0 as minimum for genome browser baseline
            data_min = data["value"].min()
            self.y_min = (
                self.min_value
                if self.min_value is not None
                else min(0, data_min)
            )
            self.y_max = (
                self.max_value
                if self.max_value is not None
                else data["value"].max()
            )

        # Ensure some range
        if self.y_min == self.y_max:
            self.y_max = self.y_min + 1

        ax.set_ylim(ymin=self.y_min, ymax=self.y_max)

        # Add labels
        if self.label.plot_title or self.label.plot_scale:
            labeller = TrackLabeller.from_config(
                self.label,
                gr,
                self.y_min,
                self.y_max,
                title=self.title or "",
            )
            labeller.plot(ax, gr)
        else:
            clean_axis(ax)
