"""Difference track for two BigWig signals."""

import matplotlib.axes
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from .base import Track, TrackLabeller
from .bigwig import BigWigTrack
from .enums import BigWigDiffMethod
from .region import GenomicRegion
from .utils import clean_axis


class BigWigDiffAesthetics(BaseModel):
    """Visual configuration for BigWigDiff tracks."""

    positive_color: str = Field(default="#d62728", description="Bar color for positive differences.")
    negative_color: str = Field(default="#1f77b4", description="Bar color for negative differences.")
    linewidth: float = Field(default=0.7, description="Edge line width for diff bars.")
    bar_alpha: float = Field(default=0.45, description="Opacity of diff bars (0-1).")
    zero_line_color: str = Field(default="#333333", description="Color of the horizontal zero baseline.")
    zero_line_width: float = Field(default=0.6, description="Line width for the zero baseline.")
    zero_line_alpha: float = Field(default=0.8, description="Opacity of the zero baseline.")

    model_config = ConfigDict(use_enum_values=True)


class BigWigDiff(Track):
    file_a: str = Field(description="Path to first BigWig input.")
    file_b: str = Field(description="Path to second BigWig input.")
    method: BigWigDiffMethod = Field(
        default=BigWigDiffMethod.SUBTRACT,
        description="Computation used to derive differential signal values.",
    )
    aesthetics: BigWigDiffAesthetics = Field(
        default_factory=BigWigDiffAesthetics,
        description="Visual styling for positive/negative diff bars and zero line.",
    )
    height: float = Field(default=1.5, description="Relative panel height for this track.")

    def _align(self, a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
        merged = pd.merge(a, b, on=["chrom", "start", "end"], suffixes=("_a", "_b"))
        return merged

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:
        track_a = BigWigTrack(data=self.file_a)
        track_b = BigWigTrack(data=self.file_b)
        a = track_a.fetch_data(gr)
        b = track_b.fetch_data(gr)

        merged = self._align(a, b)
        if merged.empty:
            return pd.DataFrame(columns=["x", "start", "end", "value"])

        values_a = merged["value_a"].to_numpy(dtype=float)
        values_b = merged["value_b"].to_numpy(dtype=float)

        if self.method == "ratio":
            value = np.divide(values_a, values_b, out=np.zeros_like(values_a), where=values_b != 0)
        elif self.method == "log2ratio":
            ratio = np.divide(values_a, values_b, out=np.ones_like(values_a), where=values_b != 0)
            value = np.log2(np.clip(ratio, 1e-12, None))
        else:
            value = values_a - values_b

        x = ((merged["start"] + merged["end"]) / 2).to_numpy(dtype=float)
        return pd.DataFrame(
            {
                "x": x,
                "start": merged["start"].to_numpy(dtype=float),
                "end": merged["end"].to_numpy(dtype=float),
                "value": value,
            }
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        data = self.fetch_data(gr)
        if data.empty:
            ax.set_xlim(gr.start, gr.end)
            ax.set_ylim(-1, 1)
            clean_axis(ax)
            return

        starts = data["start"].to_numpy(dtype=float)
        widths = (data["end"] - data["start"]).to_numpy(dtype=float)
        y = data["value"].to_numpy()

        positive = np.where(y >= 0, y, 0.0)
        negative = np.where(y < 0, y, 0.0)

        ax.bar(
            starts,
            positive,
            width=widths,
            align="edge",
            color=self.positive_color,
            edgecolor=self.positive_color,
            linewidth=self.linewidth,
            alpha=self.bar_alpha,
        )
        ax.bar(
            starts,
            negative,
            width=widths,
            align="edge",
            color=self.negative_color,
            edgecolor=self.negative_color,
            linewidth=self.linewidth,
            alpha=self.bar_alpha,
        )
        ax.axhline(
            0,
            color=self.zero_line_color,
            linewidth=self.zero_line_width,
            alpha=self.zero_line_alpha,
        )

        y_min = float(np.nanmin(y))
        y_max = float(np.nanmax(y))
        if y_min == y_max:
            y_max = y_min + 1

        TrackLabeller.from_config(
            self.label,
            gr,
            y_min,
            y_max,
            title=self.title or "",
            title_color=self.positive_color,
        ).plot(ax, gr)
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
