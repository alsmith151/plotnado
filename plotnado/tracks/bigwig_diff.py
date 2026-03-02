"""Difference track for two BigWig signals."""

from typing import Literal

import matplotlib.axes
import numpy as np
import pandas as pd
from pydantic import BaseModel

from .base import Track, TrackLabeller
from .bigwig import BigWigTrack
from .region import GenomicRegion
from .utils import clean_axis


class BigWigDiffAesthetics(BaseModel):
    """Visual configuration for BigWigDiff tracks."""

    positive_color: str = "#d62728"
    negative_color: str = "#1f77b4"
    linewidth: float = 1.0
    bar_alpha: float = 0.45
    zero_line_color: str = "#333333"
    zero_line_width: float = 0.8
    zero_line_alpha: float = 0.8


class BigWigDiff(Track):
    file_a: str
    file_b: str
    method: Literal["subtract", "ratio", "log2ratio"] = "subtract"
    aesthetics: BigWigDiffAesthetics = BigWigDiffAesthetics()
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

        TrackLabeller(
            gr=gr,
            y_min=y_min,
            y_max=y_max,
            title=self.title or self.method,
            plot_title=self.label.plot_title and bool(self.title),
            plot_scale=self.label.plot_scale,
            label_on_track=self.label.label_on_track,
            data_range_style=self.label.data_range_style,
            label_box_enabled=self.label.label_box_enabled,
            label_box_alpha=self.label.label_box_alpha,
            title_location=self.label.title_location,
            title_height=self.label.title_height,
            title_size=self.label.title_size,
            title_color=self.label.title_color,
            title_font=self.label.title_font,
            title_weight=self.label.title_weight,
            scale_location=self.label.scale_location,
            scale_height=self.label.scale_height,
            scale_precision=self.label.scale_precision,
            scale_size=self.label.scale_size,
            scale_color=self.label.scale_color,
            scale_font=self.label.scale_font,
            scale_weight=self.label.scale_weight,
        ).plot(ax, gr)
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
