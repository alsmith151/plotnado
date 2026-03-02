"""Cooler-based matrix tracks."""

from typing import Literal

import matplotlib.axes
import numpy as np
from pydantic import BaseModel

from .base import Track
from .region import GenomicRegion
from .utils import clean_axis


class CoolerAesthetics(BaseModel):
    cmap: str = "RdBu_r"
    min_value: float | None = None
    max_value: float | None = None

    @property
    def vmin(self) -> float | None:
        return self.min_value

    @vmin.setter
    def vmin(self, value: float | None) -> None:
        self.min_value = value

    @property
    def vmax(self) -> float | None:
        return self.max_value

    @vmax.setter
    def vmax(self, value: float | None) -> None:
        self.max_value = value


class CoolerTrack(Track):
    file: str
    resolution: int | None = None
    balance: bool = True
    transform: Literal["log", "log2", "log10", "none"] = "none"
    aesthetics: CoolerAesthetics = CoolerAesthetics()
    height: float = 2.5

    def _cooler_uri(self) -> str:
        if self.file.endswith(".mcool") and self.resolution:
            return f"{self.file}::resolutions/{self.resolution}"
        return self.file

    def _transform(self, matrix: np.ndarray) -> np.ndarray:
        data = np.asarray(matrix, dtype=float)
        if self.transform == "log":
            return np.log(np.clip(data, 1e-12, None))
        if self.transform == "log2":
            return np.log2(np.clip(data, 1e-12, None))
        if self.transform == "log10":
            return np.log10(np.clip(data, 1e-12, None))
        return data

    def fetch_data(self, gr: GenomicRegion) -> np.ndarray:
        try:
            import cooler
        except ImportError as exc:
            raise ImportError(
                "CoolerTrack requires optional dependency 'cooler'. Install with plotnado[cooler]."
            ) from exc

        region = f"{gr.chromosome}:{gr.start}-{gr.end}"
        matrix = cooler.Cooler(self._cooler_uri()).matrix(balance=self.balance).fetch(region)
        return self._transform(matrix)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        matrix = self.fetch_data(gr)
        if matrix.size == 0:
            clean_axis(ax)
            return

        tri = np.triu(matrix)
        ax.matshow(
            tri,
            cmap=self.cmap,
            vmin=self.min_value,
            vmax=self.max_value,
            origin="lower",
            aspect="auto",
        )
        clean_axis(ax)


class CapcruncherTrack(CoolerTrack):
    viewpoint: str | None = None
    normalisation: str | None = None


class CoolerAverage(Track):
    files: list[str]
    resolution: int | None = None
    balance: bool = True
    transform: Literal["log", "log2", "log10", "none"] = "none"
    aesthetics: CoolerAesthetics = CoolerAesthetics()
    height: float = 2.5

    def fetch_data(self, gr: GenomicRegion) -> np.ndarray:
        matrices = []
        for path in self.files:
            track = CoolerTrack(
                file=path,
                resolution=self.resolution,
                balance=self.balance,
                transform=self.transform,
            )
            matrices.append(track.fetch_data(gr))

        if not matrices:
            return np.array([])
        return np.mean(np.stack(matrices, axis=0), axis=0)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        matrix = self.fetch_data(gr)
        if matrix.size == 0:
            clean_axis(ax)
            return
        ax.matshow(
            np.triu(matrix),
            cmap=self.cmap,
            vmin=self.min_value,
            vmax=self.max_value,
            origin="lower",
            aspect="auto",
        )
        clean_axis(ax)
