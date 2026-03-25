"""QuantNado-backed tracks for multi-omics locus visualization."""

from __future__ import annotations

import importlib
from typing import Any

import matplotlib.axes
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, model_validator

from .base import Track, TrackLabeller
from .enums import TrackType
from .region import GenomicRegion
from .registry import registry
from .utils import clean_axis


def _region_string(gr: GenomicRegion) -> str:
    return f"{gr.chromosome}:{gr.start}-{gr.end}"


def _compute_if_needed(data: Any) -> Any:
    if hasattr(data, "compute"):
        return data.compute()
    return data


def _array_values(data: Any) -> np.ndarray:
    raw = getattr(data, "values", data)
    return np.asarray(raw)


def _array_positions(data: Any, default_start: int, size: int) -> np.ndarray:
    coords = getattr(data, "coords", None)
    if coords is not None and "position" in coords:
        coord = coords["position"]
        values = np.asarray(getattr(coord, "values", coord), dtype=float).reshape(-1)
        if values.size == size:
            return values
    return np.arange(size, dtype=float) + float(default_start)


def _select_sample_and_region(data: Any, sample: str, gr: GenomicRegion) -> Any:
    selected = data
    if hasattr(selected, "sel"):
        try:
            selected = selected.sel(sample=sample)
        except Exception:
            pass
        try:
            selected = selected.sel(position=slice(gr.start, gr.end))
        except Exception:
            pass
    return _compute_if_needed(selected)


def _extract_series(data: Any, sample: str, gr: GenomicRegion) -> tuple[np.ndarray, np.ndarray]:
    selected = _select_sample_and_region(data, sample, gr)
    values = _array_values(selected).astype(float).reshape(-1)
    positions = _array_positions(selected, default_start=gr.start, size=values.size)
    return positions, values


class QuantNadoCoverageAesthetics(BaseModel):
    """Aesthetics for QuantNado coverage tracks.

    Example:
        >>> QuantNadoCoverageAesthetics(color="#2171b5", fill=True)
    """
    color: str = Field(default="#2171b5", description="Primary color for coverage rendering.")
    alpha: float = Field(default=0.75, description="Opacity for coverage fill and line.")
    fill: bool = Field(default=True, description="Fill area under the coverage profile.")
    linewidth: float = Field(default=0.8, description="Line width for stepped coverage trace.")
    min_value: float | None = Field(default=None, description="Optional fixed lower y-limit.")
    max_value: float | None = Field(default=None, description="Optional fixed upper y-limit.")
    show_baseline: bool = Field(default=True, description="Draw baseline at y=0 when visible.")
    baseline_color: str = Field(default="#b8bec8", description="Baseline color.")
    baseline_alpha: float = Field(default=0.85, description="Baseline opacity.")
    baseline_linewidth: float = Field(default=0.6, description="Baseline line width.")


class QuantNadoStrandedCoverageAesthetics(BaseModel):
    """Aesthetics for QuantNado stranded coverage tracks.

    Example:
        >>> QuantNadoStrandedCoverageAesthetics(color="#1f78b4", reverse_color="#d62728")
    """
    color: str = Field(default="#1f78b4", description="Color for forward strand signal.")
    reverse_color: str | None = Field(
        default=None,
        description="Optional color override for reverse strand signal; defaults to `color`.",
    )
    alpha: float = Field(default=0.7, description="Opacity for stranded coverage traces.")
    fill: bool = Field(default=True, description="Fill areas for forward and reverse signals.")
    linewidth: float = Field(default=0.8, description="Line width for stepped traces.")
    min_value: float | None = Field(default=None, description="Optional fixed lower y-limit.")
    max_value: float | None = Field(default=None, description="Optional fixed upper y-limit.")
    show_baseline: bool = Field(default=True, description="Draw baseline at y=0.")
    baseline_color: str = Field(default="#444444", description="Baseline color.")
    baseline_alpha: float = Field(default=0.8, description="Baseline opacity.")
    baseline_linewidth: float = Field(default=0.6, description="Baseline line width.")


class QuantNadoMethylationAesthetics(BaseModel):
    """Aesthetics for QuantNado methylation tracks.

    Example:
        >>> QuantNadoMethylationAesthetics(color="#2a9d8f", point_size=12.0)
    """
    color: str = Field(default="#2a9d8f", description="Scatter color for methylation points.")
    alpha: float = Field(default=0.75, description="Opacity for methylation points.")
    point_size: float = Field(default=10.0, description="Marker area for methylation points.")
    min_value: float | None = Field(default=0.0, description="Lower y-limit (default 0).")
    max_value: float | None = Field(default=100.0, description="Upper y-limit (default 100).")


class QuantNadoVariantAesthetics(BaseModel):
    """Aesthetics for QuantNado variant tracks.

    Example:
        >>> QuantNadoVariantAesthetics(het_color="#1f77b4", hom_alt_color="#d62728")
    """
    het_color: str = Field(default="#1f77b4", description="Color for heterozygous variants.")
    hom_alt_color: str = Field(default="#d62728", description="Color for homozygous-alt variants.")
    alpha: float = Field(default=0.8, description="Opacity for lollipop stems/markers.")
    linewidth: float = Field(default=0.9, description="Lollipop stem line width.")
    marker_size: float = Field(default=24.0, description="Lollipop marker area.")
    min_value: float | None = Field(default=0.0, description="Lower y-limit (default 0).")
    max_value: float | None = Field(default=1.08, description="Upper y-limit (default 1.08).")
    show_baseline: bool = Field(default=True, description="Draw baseline at y=0.")
    baseline_color: str = Field(default="#444444", description="Baseline color.")
    baseline_alpha: float = Field(default=0.75, description="Baseline opacity.")
    baseline_linewidth: float = Field(default=0.6, description="Baseline line width.")


class _QuantNadoSourceMixin:
    quantnado: Any | None = Field(
        default=None,
        description="Runtime QuantNado instance used to fetch track data.",
    )
    dataset_path: str | None = Field(
        default=None,
        description="Path to a QuantNado dataset opened lazily when needed.",
    )

    _qn_instance: Any | None = PrivateAttr(default=None)

    def _resolve_quantnado(self) -> Any | None:
        if self.quantnado is not None:
            return self.quantnado
        if self._qn_instance is not None:
            return self._qn_instance
        if self.dataset_path is None:
            return None

        try:
            module = importlib.import_module("quantnado")
        except ImportError as exc:
            raise ImportError(
                "dataset_path requires optional dependency 'quantnado'. "
                "Install with plotnado[quantnado]."
            ) from exc

        opener = getattr(module, "open", None)
        if callable(opener):
            self._qn_instance = opener(self.dataset_path)
            return self._qn_instance

        quantnado_cls = getattr(module, "QuantNado", None)
        if quantnado_cls is None or not hasattr(quantnado_cls, "open"):
            raise ImportError("Unable to locate quantnado.open or QuantNado.open")

        self._qn_instance = quantnado_cls.open(self.dataset_path)
        return self._qn_instance

    @staticmethod
    def _label_or_clean(
        ax: matplotlib.axes.Axes,
        track: Track,
        gr: GenomicRegion,
        y_min: float,
        y_max: float,
        title_color: str | None = None,
    ) -> None:
        if track.label.plot_title or track.label.plot_scale:
            TrackLabeller.from_config(
                track.label,
                gr,
                y_min,
                y_max,
                title=track.title or "",
                title_color=title_color,
            ).plot(ax, gr)
        else:
            clean_axis(ax)


@registry.register(TrackType.QUANTNADO_COVERAGE)
class QuantNadoCoverageTrack(_QuantNadoSourceMixin, Track):
    """Track for plotting QuantNado per-base or binned coverage.

    Example:
        >>> QuantNadoCoverageTrack(sample="tumor", dataset_path="study.qn")
    """
    sample: str = Field(description="Sample name to render from QuantNado data.")
    scaling_factor: float = Field(
        default=1.0,
        description="User-defined multiplicative factor applied to the extracted coverage signal.",
    )
    normalise: str | None = Field(
        default=None,
        description="Optional normalization mode passed to QuantNado extract_region, e.g. 'cpm' or 'rpkm'.",
    )
    normalize: str | None = Field(
        default=None,
        description="American-English alias for `normalise`.",
    )
    library_sizes: pd.Series | dict | None = Field(
        default=None,
        description="Optional library sizes forwarded to QuantNado normalization.",
    )
    coverage_data: Any | None = Field(
        default=None,
        description="Optional precomputed xarray-like coverage data with dims (sample, position).",
    )
    aesthetics: QuantNadoCoverageAesthetics = Field(
        default_factory=QuantNadoCoverageAesthetics,
        description="Visual styling options for coverage rendering.",
    )

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def _validate_source(self) -> "QuantNadoCoverageTrack":
        if self.quantnado is None and self.dataset_path is None and self.coverage_data is None:
            raise ValueError(
                "QuantNadoCoverageTrack requires one source: quantnado, dataset_path, or coverage_data."
            )
        return self

    def fetch_data(self, gr: GenomicRegion) -> dict[str, np.ndarray]:
        if self.coverage_data is not None:
            positions, values = _extract_series(self.coverage_data, self.sample, gr)
            values = values * float(self.scaling_factor)
            return {"position": positions, "value": values}

        qn = self._resolve_quantnado()
        if qn is None:
            raise RuntimeError("No QuantNado source available for coverage track")
        coverage_store = getattr(qn, "coverage", None)
        if coverage_store is None:
            raise RuntimeError("QuantNado source has no coverage store for coverage track")
        data = coverage_store.extract_region(
            region=_region_string(gr),
            samples=[self.sample],
            as_xarray=False,
            normalise=self.normalise,
            normalize=self.normalize,
            library_sizes=self.library_sizes,
        )
        positions, values = _extract_series(data, self.sample, gr)
        values = values * float(self.scaling_factor)
        return {"position": positions, "value": values}

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        fetched = self.fetch_data(gr)
        x = fetched["position"]
        y = fetched["value"]

        if x.size and y.size:
            if self.fill:
                ax.fill_between(x, y, step="post", alpha=self.alpha, color=self.color)
            ax.plot(
                x,
                y,
                drawstyle="steps-post",
                color=self.color,
                linewidth=self.linewidth,
                alpha=self.alpha,
            )

        if y.size:
            finite = y[np.isfinite(y)]
            if finite.size:
                y_min = float(np.nanmin(finite))
                y_max = float(np.nanmax(finite))
            else:
                y_min, y_max = 0.0, 1.0
        else:
            y_min, y_max = 0.0, 1.0

        if self.min_value is not None:
            y_min = float(self.min_value)
        else:
            y_min = min(0.0, y_min)
        if self.max_value is not None:
            y_max = float(self.max_value)
        if y_min == y_max:
            y_max = y_min + 1.0

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
        if self.show_baseline and y_min <= 0 <= y_max:
            ax.axhline(
                0,
                color=self.baseline_color,
                alpha=self.baseline_alpha,
                linewidth=self.baseline_linewidth,
                zorder=0,
            )

        self._label_or_clean(ax, self, gr, y_min, y_max, title_color=self.color)


@registry.register(TrackType.QUANTNADO_STRANDED_COVERAGE)
class QuantNadoStrandedCoverageTrack(_QuantNadoSourceMixin, Track):
    """Track for plotting forward and reverse QuantNado coverage together.

    Example:
        >>> QuantNadoStrandedCoverageTrack(sample="tumor", dataset_path="study.qn")
    """
    sample: str = Field(description="Sample name to render from QuantNado data.")
    scaling_factor: float = Field(
        default=1.0,
        description="User-defined multiplicative factor applied to forward and reverse coverage signals.",
    )
    normalise: str | None = Field(
        default=None,
        description="Optional normalization mode passed to QuantNado coverage store extraction.",
    )
    normalize: str | None = Field(
        default=None,
        description="American-English alias for `normalise`.",
    )
    library_sizes: pd.Series | dict | None = Field(
        default=None,
        description="Optional library sizes forwarded to QuantNado normalization.",
    )
    coverage_fwd_data: Any | None = Field(
        default=None,
        description="Optional forward-strand xarray-like coverage data with dims (sample, position).",
    )
    coverage_rev_data: Any | None = Field(
        default=None,
        description="Optional reverse-strand xarray-like coverage data with dims (sample, position).",
    )
    aesthetics: QuantNadoStrandedCoverageAesthetics = Field(
        default_factory=QuantNadoStrandedCoverageAesthetics,
        description="Visual styling options for stranded coverage rendering.",
    )

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def _validate_source(self) -> "QuantNadoStrandedCoverageTrack":
        has_object_source = self.quantnado is not None or self.dataset_path is not None
        has_array_source = self.coverage_fwd_data is not None and self.coverage_rev_data is not None
        if not has_object_source and not has_array_source:
            raise ValueError(
                "QuantNadoStrandedCoverageTrack requires quantnado/dataset_path "
                "or both coverage_fwd_data and coverage_rev_data."
            )
        return self

    def fetch_data(self, gr: GenomicRegion) -> dict[str, np.ndarray]:
        if self.coverage_fwd_data is not None and self.coverage_rev_data is not None:
            pos_fwd, fwd = _extract_series(self.coverage_fwd_data, self.sample, gr)
            pos_rev, rev = _extract_series(self.coverage_rev_data, self.sample, gr)
            factor = float(self.scaling_factor)
            fwd = fwd * factor
            rev = rev * factor
            if pos_fwd.size and pos_rev.size and np.array_equal(pos_fwd, pos_rev):
                pos = pos_fwd
            elif pos_fwd.size:
                pos = pos_fwd
            else:
                pos = pos_rev
            return {"position": pos, "forward": fwd, "reverse": rev}

        qn = self._resolve_quantnado()
        if qn is None:
            raise RuntimeError("No QuantNado source available for stranded coverage track")
        coverage_store = getattr(qn, "coverage", None)
        if coverage_store is None:
            raise RuntimeError("QuantNado source has no coverage store for stranded coverage track")

        region = _region_string(gr)
        fwd_data = coverage_store.extract_region(
            region=region,
            samples=[self.sample],
            as_xarray=False,
            strand="+",
            normalise=self.normalise,
            normalize=self.normalize,
            library_sizes=self.library_sizes,
        )
        rev_data = coverage_store.extract_region(
            region=region,
            samples=[self.sample],
            as_xarray=False,
            strand="-",
            normalise=self.normalise,
            normalize=self.normalize,
            library_sizes=self.library_sizes,
        )
        pos_fwd, fwd = _extract_series(fwd_data, self.sample, gr)
        pos_rev, rev = _extract_series(rev_data, self.sample, gr)
        factor = float(self.scaling_factor)
        fwd = fwd * factor
        rev = rev * factor
        if pos_fwd.size and pos_rev.size and np.array_equal(pos_fwd, pos_rev):
            pos = pos_fwd
        elif pos_fwd.size:
            pos = pos_fwd
        else:
            pos = pos_rev
        return {"position": pos, "forward": fwd, "reverse": rev}

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        fetched = self.fetch_data(gr)
        x = fetched["position"]
        fwd = fetched["forward"]
        rev = fetched["reverse"]

        reverse_color = self.reverse_color or self.color

        if x.size:
            if self.fill:
                ax.fill_between(x, fwd, step="post", alpha=self.alpha, color=self.color)
                ax.fill_between(x, -rev, step="post", alpha=self.alpha, color=reverse_color)
            ax.plot(
                x,
                fwd,
                drawstyle="steps-post",
                color=self.color,
                linewidth=self.linewidth,
                alpha=self.alpha,
            )
            ax.plot(
                x,
                -rev,
                drawstyle="steps-post",
                color=reverse_color,
                linewidth=self.linewidth,
                alpha=self.alpha,
            )

        max_fwd = float(np.nanmax(np.abs(fwd))) if fwd.size else 0.0
        max_rev = float(np.nanmax(np.abs(rev))) if rev.size else 0.0
        y_extent = max(max_fwd, max_rev, 1.0)
        y_min = -y_extent
        y_max = y_extent
        if self.min_value is not None:
            y_min = float(self.min_value)
        if self.max_value is not None:
            y_max = float(self.max_value)
        if y_min == y_max:
            y_max = y_min + 1.0

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
        if self.show_baseline and y_min <= 0 <= y_max:
            ax.axhline(
                0,
                color=self.baseline_color,
                alpha=self.baseline_alpha,
                linewidth=self.baseline_linewidth,
                zorder=0,
            )

        self._label_or_clean(ax, self, gr, y_min, y_max, title_color=self.color)


@registry.register(TrackType.QUANTNADO_METHYLATION)
class QuantNadoMethylationTrack(_QuantNadoSourceMixin, Track):
    """Track for plotting QuantNado methylation measurements.

    Example:
        >>> QuantNadoMethylationTrack(sample="tumor", dataset_path="study.qn")
    """
    sample: str = Field(description="Sample name to render from QuantNado data.")
    methylation_variable: str = Field(
        default="methylation_pct",
        description="Variable fetched from QuantNado methylation store.",
    )
    methylation_data: Any | None = Field(
        default=None,
        description="Optional precomputed methylation xarray-like data with dims (sample, position).",
    )
    aesthetics: QuantNadoMethylationAesthetics = Field(
        default_factory=QuantNadoMethylationAesthetics,
        description="Visual styling options for methylation rendering.",
    )

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def _validate_source(self) -> "QuantNadoMethylationTrack":
        if self.quantnado is None and self.dataset_path is None and self.methylation_data is None:
            raise ValueError(
                "QuantNadoMethylationTrack requires one source: quantnado, dataset_path, or methylation_data."
            )
        return self

    def fetch_data(self, gr: GenomicRegion) -> dict[str, np.ndarray]:
        if self.methylation_data is not None:
            positions, values = _extract_series(self.methylation_data, self.sample, gr)
            return {"position": positions, "value": values}

        qn = self._resolve_quantnado()
        if qn is None:
            raise RuntimeError("No QuantNado source available for methylation track")
        methylation_store = getattr(qn, "methylation", None)
        if methylation_store is None:
            raise RuntimeError("QuantNado source has no methylation store")
        data = methylation_store.extract_region(
            region=_region_string(gr),
            variable=self.methylation_variable,
            samples=[self.sample],
        )
        positions, values = _extract_series(data, self.sample, gr)
        return {"position": positions, "value": values}

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        fetched = self.fetch_data(gr)
        x = fetched["position"]
        y = fetched["value"]

        if x.size and y.size:
            ax.scatter(
                x,
                y,
                s=self.point_size,
                color=self.color,
                alpha=self.alpha,
                linewidths=0,
                zorder=2,
            )

        y_min = float(self.min_value) if self.min_value is not None else 0.0
        y_max = float(self.max_value) if self.max_value is not None else 100.0
        if y_min == y_max:
            y_max = y_min + 1.0

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
        self._label_or_clean(ax, self, gr, y_min, y_max, title_color=self.color)


@registry.register(TrackType.QUANTNADO_VARIANT)
class QuantNadoVariantTrack(_QuantNadoSourceMixin, Track):
    """Track for plotting QuantNado variant allele fractions.

    Example:
        >>> QuantNadoVariantTrack(sample="tumor", dataset_path="study.qn")
    """
    sample: str = Field(description="Sample name to render from QuantNado data.")
    allele_depth_ref_variable: str = Field(
        default="allele_depth_ref",
        description="Reference-allele depth variable name in QuantNado variants store.",
    )
    allele_depth_alt_variable: str = Field(
        default="allele_depth_alt",
        description="Alternate-allele depth variable name in QuantNado variants store.",
    )
    genotype_variable: str = Field(
        default="genotype",
        description="Genotype variable name in QuantNado variants store.",
    )
    fetch_genotype: bool = Field(
        default=True,
        description="Fetch genotype data from QuantNado variants store when available.",
    )
    allele_depth_ref_data: Any | None = Field(
        default=None,
        description="Optional precomputed ref-depth xarray-like data with dims (sample, position).",
    )
    allele_depth_alt_data: Any | None = Field(
        default=None,
        description="Optional precomputed alt-depth xarray-like data with dims (sample, position).",
    )
    genotype_data: Any | None = Field(
        default=None,
        description="Optional precomputed genotype xarray-like data with dims (sample, position).",
    )
    aesthetics: QuantNadoVariantAesthetics = Field(
        default_factory=QuantNadoVariantAesthetics,
        description="Visual styling options for variant lollipop rendering.",
    )

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def _validate_source(self) -> "QuantNadoVariantTrack":
        has_object_source = self.quantnado is not None or self.dataset_path is not None
        has_array_source = self.allele_depth_ref_data is not None and self.allele_depth_alt_data is not None
        if not has_object_source and not has_array_source:
            raise ValueError(
                "QuantNadoVariantTrack requires quantnado/dataset_path "
                "or both allele_depth_ref_data and allele_depth_alt_data."
            )
        return self

    @staticmethod
    def _classify_genotype(af: np.ndarray, provided: np.ndarray | None) -> np.ndarray:
        if provided is not None:
            return provided.astype(int)
        gt = np.full(af.shape, -1, dtype=int)
        gt[(af > 0.2) & (af < 0.8)] = 1
        gt[af >= 0.8] = 2
        return gt

    def fetch_data(self, gr: GenomicRegion) -> dict[str, np.ndarray]:
        if self.allele_depth_ref_data is not None and self.allele_depth_alt_data is not None:
            positions, ref_depth = _extract_series(self.allele_depth_ref_data, self.sample, gr)
            _, alt_depth = _extract_series(self.allele_depth_alt_data, self.sample, gr)
            genotype = None
            if self.genotype_data is not None:
                _, genotype_raw = _extract_series(self.genotype_data, self.sample, gr)
                genotype = genotype_raw.astype(int)
            return {
                "position": positions,
                "ref_depth": ref_depth,
                "alt_depth": alt_depth,
                "genotype": genotype,
            }

        qn = self._resolve_quantnado()
        if qn is None:
            raise RuntimeError("No QuantNado source available for variant track")
        variants_store = getattr(qn, "variants", None)
        if variants_store is None:
            raise RuntimeError("QuantNado source has no variants store")

        region = _region_string(gr)
        ref_data = variants_store.extract_region(
            region=region,
            variable=self.allele_depth_ref_variable,
            samples=[self.sample],
        )
        alt_data = variants_store.extract_region(
            region=region,
            variable=self.allele_depth_alt_variable,
            samples=[self.sample],
        )
        positions, ref_depth = _extract_series(ref_data, self.sample, gr)
        _, alt_depth = _extract_series(alt_data, self.sample, gr)

        genotype = None
        if self.fetch_genotype:
            try:
                genotype_data = variants_store.extract_region(
                    region=region,
                    variable=self.genotype_variable,
                    samples=[self.sample],
                )
                _, genotype_raw = _extract_series(genotype_data, self.sample, gr)
                genotype = genotype_raw.astype(int)
            except Exception:
                genotype = None

        return {
            "position": positions,
            "ref_depth": ref_depth,
            "alt_depth": alt_depth,
            "genotype": genotype,
        }

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        fetched = self.fetch_data(gr)
        x = fetched["position"]
        ref_depth = fetched["ref_depth"]
        alt_depth = fetched["alt_depth"]
        genotype_raw = fetched["genotype"]

        total = ref_depth + alt_depth
        with np.errstate(invalid="ignore", divide="ignore"):
            af = np.where(total > 0, alt_depth / total, np.nan)
        genotype = self._classify_genotype(af, genotype_raw)

        for genotype_value, color in ((1, self.het_color), (2, self.hom_alt_color)):
            mask = genotype == genotype_value
            if not np.any(mask):
                continue
            ax.vlines(
                x[mask],
                0,
                af[mask],
                color=color,
                linewidth=self.linewidth,
                alpha=self.alpha,
            )
            ax.scatter(
                x[mask],
                af[mask],
                color=color,
                s=self.marker_size,
                alpha=self.alpha,
                zorder=3,
            )

        y_min = float(self.min_value) if self.min_value is not None else 0.0
        y_max = float(self.max_value) if self.max_value is not None else 1.08
        if y_min == y_max:
            y_max = y_min + 1.0

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(y_min, y_max)
        if self.show_baseline and y_min <= 0 <= y_max:
            ax.axhline(
                0,
                color=self.baseline_color,
                alpha=self.baseline_alpha,
                linewidth=self.baseline_linewidth,
                zorder=0,
            )

        self._label_or_clean(ax, self, gr, y_min, y_max, title_color=self.het_color)
