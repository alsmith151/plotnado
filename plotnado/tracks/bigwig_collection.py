"""Collection track for rendering multiple BigWig files."""

from pathlib import Path

import matplotlib.axes
import pandas as pd
from pydantic import BaseModel, ConfigDict

from .base import LabelConfig, Track, TrackLabeller
from .bigwig import BigWigTrack, BigwigAesthetics
from .region import GenomicRegion
from .scaling import Autoscaler
from .utils import clean_axis
from .enums import CollectionStyle


class BigWigCollectionAesthetics(BaseModel):
    colors: list[str] | None = None
    labels: list[str] | None = None
    alpha: float = 0.6
    style: CollectionStyle = CollectionStyle.OVERLAY

    model_config = ConfigDict(use_enum_values=True)


class BigWigCollection(Track):
    files: list[str]
    aesthetics: BigWigCollectionAesthetics = BigWigCollectionAesthetics()
    height: float = 2.0

    def _build_tracks(self) -> list[BigWigTrack]:
        defaults = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
        tracks: list[BigWigTrack] = []
        for index, file_path in enumerate(self.files):
            color = (
                self.colors[index]
                if self.colors and index < len(self.colors)
                else defaults[index % len(defaults)]
            )
            label = (
                self.labels[index]
                if self.labels and index < len(self.labels)
                else Path(file_path).stem
            )
            tracks.append(
                BigWigTrack(
                    data=file_path,
                    title=label,
                    aesthetics=BigwigAesthetics(
                        style="fill",
                        color=color,
                        alpha=self.alpha,
                    ),
                    label=LabelConfig(plot_title=False, plot_scale=False),
                )
            )
        return tracks

    def fetch_data(self, gr: GenomicRegion) -> list[pd.DataFrame]:
        return [track.fetch_data(gr) for track in self._build_tracks()]

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        tracks = self._build_tracks()
        scaler = Autoscaler(tracks=tracks, gr=gr)
        y_min = scaler.min_value
        y_max = scaler.max_value
        if y_min == y_max:
            y_max = y_min + 1

        if self.style == "overlay":
            for track in tracks:
                track.min_value = y_min
                track.max_value = y_max
                track.plot(ax, gr)
        else:
            offset = 0.0
            max_step = y_max - y_min
            for track in tracks:
                data = track.fetch_data(gr)
                if data.empty:
                    continue
                x = (data["start"] + data["end"]) / 2
                y = data["value"] + offset
                ax.plot(x, y, color=track.color, alpha=track.alpha)
                offset += max_step
            y_max = max_step * max(1, len(tracks))
            y_min = 0

        if self.title:
            TrackLabeller(
                gr=gr,
                y_min=y_min,
                y_max=y_max,
                title=self.title,
                plot_title=True,
                plot_scale=True,
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
        else:
            clean_axis(ax)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
