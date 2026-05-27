"""
BigWig track for signal visualization.
"""

from pathlib import Path
from typing import Literal

import matplotlib.axes
import matplotlib.patches
import numpy as np
import pandas as pd
import pybigtools
from pydantic import Field

from .region import GenomicRegion
from .base import Track, TrackLabeller
from .utils import clean_axis
from .schemas import BedgraphDataFrame
from .enums import PlotStyle, TrackType
from .aesthetics import BaseAesthetics
from .registry import registry


class BigwigAesthetics(BaseAesthetics):
    """
    Aesthetics configuration for BigWig tracks.

    Inherits color, alpha, and linewidth from BaseAesthetics.

    Attributes:
        style: Plot style ('std', 'line', 'scatter', 'heatmap')
    """

    style: PlotStyle = Field(
        default=PlotStyle.STD,
        description="Style for rendering the signal track, for options see PlotStyle enum.",
    )
    fill: bool = Field(default=True, description="Fill the area under the signal curve when supported.")
    scatter_point_size: float = Field(
        default=1.0,
        description="Marker area for scatter style rendering.",
    )
    show_baseline: bool = Field(default=True, description="Draw a subtle baseline at y=0.")
    baseline_color: str = Field(default="#b8bec8", description="Color of the y=0 signal baseline.")
    baseline_alpha: float = Field(default=0.85, description="Opacity of the y=0 signal baseline.")
    baseline_linewidth: float = Field(default=0.6, description="Line width of the y=0 signal baseline.")
    smoothing_window: int = Field(
        default=1,
        ge=1,
        description="Rolling window size in bins for optional signal smoothing.",
    )
    smoothing_method: Literal["mean", "median"] = Field(
        default="mean",
        description="Aggregation used for smoothing when smoothing_window > 1.",
    )
    smoothing_center: bool = Field(
        default=True,
        description="Center the smoothing window around each bin when enabled.",
    )

    min_value: float | None = Field(
        default=None,
        description="Optional fixed lower y-limit; auto-derived when omitted.",
    )
    max_value: float | None = Field(
        default=None,
        description="Optional fixed upper y-limit; auto-derived when omitted.",
    )

    n_bins: int | None = Field(
        default=None,
        ge=1,
        description="Divide the plotted region into this many equal bins. Overrides bin_size.",
    )
    bin_size: int | None = Field(
        default=None,
        ge=1,
        description="Bin width in base pairs. Ignored when n_bins is set.",
    )


@registry.register(TrackType.BIGWIG, aliases=["bw", "signal", "bedgraph"])
class BigWigTrack(Track):
    """
    Track for displaying BigWig signal data.

    Attributes:
        data: Path to BigWig file or DataFrame with bedgraph-like data
        aesthetics: Visual styling configuration
    """

    data: Path | pd.DataFrame | str | None = Field(
        default=None,
        description="BigWig file path or bedgraph-like DataFrame data source.",
    )
    aesthetics: BigwigAesthetics = Field(
        default_factory=BigwigAesthetics,
        description="Visual styling options for this BigWig track.",
    )

    y_min: float | None = Field(
        default=None,
        description="Computed lower y-limit for the last plotted region.",
    )
    y_max: float | None = Field(
        default=None,
        description="Computed upper y-limit for the last plotted region.",
    )

    def _plot_stairs(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, values: pd.DataFrame
    ) -> None:
        """Plot as stairs."""
        if values.empty:
            return

        ordered = values.sort_values(["start", "end"]).reset_index(drop=True)
        starts = ordered["start"].to_numpy(dtype=float)
        ends = ordered["end"].to_numpy(dtype=float)
        signal = ordered["value"].to_numpy(dtype=float)

        # Preserve true BigWig bin geometry instead of using evenly-spaced edges.
        edges = np.concatenate(([starts[0]], ends))
        ax.stairs(
            values=signal,
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

    def _resolve_n_bins(self, gr: GenomicRegion) -> int | None:
        if self.n_bins is not None:
            return int(self.n_bins)
        if self.bin_size is not None:
            return max(1, int((gr.end - gr.start) / self.bin_size))
        return None

    def _rebin(self, data: pd.DataFrame, gr: GenomicRegion, n_bins: int) -> pd.DataFrame:
        """Weighted-average rebin of bedgraph data into n_bins equal-width bins over the region."""
        if data.empty:
            return data

        edges = np.linspace(gr.start, gr.end, n_bins + 1)
        bin_starts = edges[:-1]
        bin_ends = edges[1:]

        src_s = data["start"].to_numpy(dtype=float)
        src_e = data["end"].to_numpy(dtype=float)
        src_v = data["value"].to_numpy(dtype=float)
        valid = ~np.isnan(src_v)
        src_s, src_e, src_v = src_s[valid], src_e[valid], src_v[valid]

        bin_values = np.full(n_bins, np.nan)
        for i in range(n_bins):
            overlap = np.maximum(0.0, np.minimum(src_e, bin_ends[i]) - np.maximum(src_s, bin_starts[i]))
            w = overlap.sum()
            if w > 0:
                bin_values[i] = (src_v * overlap).sum() / w

        keep = ~np.isnan(bin_values)
        return pd.DataFrame({
            "chrom": gr.chromosome,
            "start": bin_starts[keep].astype(int),
            "end": bin_ends[keep].astype(int),
            "value": bin_values[keep],
        })

    def fetch_data(self, gr: GenomicRegion) -> BedgraphDataFrame:
        """Fetch data for the given genomic region."""
        if isinstance(self.data, pd.DataFrame):
            df = self._fetch_from_df(gr)
        else:
            df = self._fetch_from_disk(gr)

        n = self._resolve_n_bins(gr)
        if n is not None:
            df = self._rebin(df, gr, n)

        return BedgraphDataFrame(df)

    def _apply_smoothing(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply optional rolling smoothing to signal values."""
        if data.empty or self.smoothing_window <= 1:
            return data

        smoothed = data.sort_values(["start", "end"]).copy()
        values = pd.to_numeric(smoothed["value"], errors="coerce")
        rolling = values.rolling(
            window=int(self.smoothing_window),
            min_periods=1,
            center=bool(self.smoothing_center),
        )
        if self.smoothing_method == "median":
            smoothed["value"] = rolling.median()
        else:
            smoothed["value"] = rolling.mean()
        return smoothed

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
        facecolor = self.color if self.fill else "none"
        edgecolor = self.color if not self.fill and self.linewidth > 0 else "none"
        linewidth = self.linewidth if not self.fill else 0.0
        ax.bar(
            starts,
            values,
            width=widths,
            align="edge",
            bottom=0,
            color=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=self.alpha,
        )
    
    def _plot_heatmap(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, data: pd.DataFrame
    ) -> None:
        """Plot as heatmap."""
        import matplotlib.pyplot as plt

        if data.empty:
            return

        x = (data["start"] + data["end"]) / 2
        y = np.zeros_like(x)
        values = data["value"].to_numpy(dtype=float)

        # Create a colormap based on the values
        norm = plt.Normalize(vmin=values.min(), vmax=values.max())
        cmap = plt.get_cmap(self.color)
        colors = cmap(norm(values))

        ax.imshow(
            colors[np.newaxis, :, :],
            aspect="auto",
            extent=(gr.start, gr.end, self.y_min, self.y_max),
            alpha=self.alpha,
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the BigWig track."""
        data = self._apply_smoothing(self.fetch_data(gr))

        # Plot based on style
        match self.style:
            case PlotStyle.STD:
                self._plot_stairs(ax, gr, data)
            case PlotStyle.FILL:
                self._plot_fill(ax, gr, data)
            case PlotStyle.LINE:
                self._plot_line(ax, gr, data)
            case PlotStyle.SCATTER:
                self._plot_scatter(ax, gr, data)
            case PlotStyle.FRAGMENT:
                self._plot_fragment(ax, gr, data)
            case PlotStyle.HEATMAP:
                self._heatmap(ax, gr, data)
            case _:
                raise ValueError(f"Unsupported plot style: {self.style}")
        
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
        if self.show_baseline and self.y_min <= 0 <= self.y_max:
            ax.axhline(
                0,
                color=self.baseline_color,
                alpha=self.baseline_alpha,
                linewidth=self.baseline_linewidth,
                zorder=0,
            )

        # Add labels
        if self.label.plot_title or self.label.plot_scale:
            labeller = TrackLabeller.from_config(
                self.label,
                gr,
                self.y_min,
                self.y_max,
                title=self.title or "",
                title_color=self.color,
            )
            labeller.plot(ax, gr)
        else:
            clean_axis(ax)
