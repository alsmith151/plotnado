"""Difference track for two BigWig signals."""

from typing import Literal

import matplotlib.axes
import numpy as np
import pandas as pd

from .base import Track, TrackLabeller
from .bigwig import BigWigTrack
from .region import GenomicRegion
from .utils import clean_axis


class BigWigDiff(Track):
    file_a: str
    file_b: str
    method: Literal["subtract", "ratio", "log2ratio"] = "subtract"
    positive_color: str = "#d62728"
    negative_color: str = "#1f77b4"
    linewidth: float = 1.0
    height: float = 1.5

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
            alpha=0.45,
        )
        ax.bar(
            starts,
            negative,
            width=widths,
            align="edge",
            color=self.negative_color,
            edgecolor=self.negative_color,
            linewidth=self.linewidth,
            alpha=0.45,
        )
        ax.axhline(0, color="#333333", linewidth=0.8, alpha=0.8)

        y_min = float(np.nanmin(y))
        y_max = float(np.nanmax(y))
        if y_min == y_max:
            y_max = y_min + 1

        TrackLabeller(
            gr=gr,
            y_min=y_min,
            y_max=y_max,
            title=self.title or self.method,
            plot_title=bool(self.title),
            plot_scale=True,
            label_on_track=self.label_on_track,
            data_range_style=self.data_range_style,
            label_box_enabled=self.label_box_enabled,
            label_box_alpha=self.label_box_alpha,
        ).plot(ax, gr)
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
