"""Cooler-based matrix tracks."""

import matplotlib.axes
import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from .base import Track
from .enums import CoolerTransform, TrackType
from .region import GenomicRegion
from .registry import registry
from .utils import clean_axis


class CoolerAesthetics(BaseModel):
    """Visual options for cooler-derived matrix tracks.

    Example:
        >>> CoolerAesthetics(cmap="magma", min_value=0.0, max_value=5.0)
    """
    cmap: str = Field(default="RdBu_r", description="Matplotlib colormap used to render matrix intensity.")
    min_value: float | None = Field(
        default=None,
        description="Optional fixed lower bound for colormap normalization.",
    )
    max_value: float | None = Field(
        default=None,
        description="Optional fixed upper bound for colormap normalization.",
    )

    model_config = ConfigDict(use_enum_values=True)

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


@registry.register(TrackType.COOLER)
class CoolerTrack(Track):
    """Track for rendering a cooler or mcool contact matrix.

    Example:
        >>> CoolerTrack(file="contacts.mcool", resolution=10000)
    """
    file: str = Field(description="Path to cooler/mcool matrix file.")
    resolution: int | None = Field(
        default=None,
        description="Resolution bin size for .mcool files.",
    )
    balance: bool = Field(default=True, description="Use balanced matrix values when available.")
    transform: CoolerTransform = Field(
        default=CoolerTransform.NONE,
        description="Transform applied to matrix values before plotting.",
    )
    aesthetics: CoolerAesthetics = Field(
        default_factory=CoolerAesthetics,
        description="Colormap and value-range styling options.",
    )
    height: float = Field(default=2.5, description="Relative panel height for this track.")

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


@registry.register(TrackType.CAPCRUNCHER)
class CapcruncherTrack(CoolerTrack):
    """Capture-centric cooler track with viewpoint metadata.

    Example:
        >>> CapcruncherTrack(file="capture.mcool", viewpoint="MYC")
    """
    viewpoint: str | None = Field(default=None, description="Optional viewpoint identifier for capture-centric views.")
    normalisation: str | None = Field(default=None, description="Optional normalization mode label.")


@registry.register(TrackType.COOLER_AVERAGE)
class CoolerAverage(Track):
    """Track for averaging multiple cooler matrices before plotting.

    Example:
        >>> CoolerAverage(files=["a.cool", "b.cool"], resolution=10000)
    """
    files: list[str] = Field(description="List of cooler/mcool files to average.")
    resolution: int | None = Field(
        default=None,
        description="Resolution bin size for .mcool input files.",
    )
    balance: bool = Field(default=True, description="Use balanced matrices where available.")
    transform: CoolerTransform = Field(
        default=CoolerTransform.NONE,
        description="Transform applied before averaging and plotting.",
    )
    aesthetics: CoolerAesthetics = Field(
        default_factory=CoolerAesthetics,
        description="Colormap and value-range styling options.",
    )
    height: float = Field(default=2.5, description="Relative panel height for this track.")

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
