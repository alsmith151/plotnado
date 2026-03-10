"""Tests for QuantNado-backed plotnado track adapters."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

xr = pytest.importorskip("xarray")

from plotnado.tracks import (
    GenomicRegion,
    QuantNadoCoverageTrack,
    QuantNadoMethylationTrack,
    QuantNadoStrandedCoverageTrack,
    QuantNadoVariantTrack,
)


def _make_dataarray(values: list[float], sample: str = "s1", start: int = 100) -> xr.DataArray:
    arr = np.asarray(values, dtype=float).reshape(1, -1)
    positions = np.arange(start, start + arr.shape[1])
    return xr.DataArray(
        arr,
        dims=("sample", "position"),
        coords={"sample": [sample], "position": positions},
    )


class _RecordingStore:
    def __init__(self, responses: dict[str, xr.DataArray] | xr.DataArray):
        self._responses = responses
        self.calls: list[dict] = []

    def extract_region(self, region: str, samples: list[str] | None = None, **kwargs):
        payload = {"region": region, "samples": samples, **kwargs}
        self.calls.append(payload)
        if isinstance(self._responses, dict):
            variable = kwargs.get("variable")
            if variable in self._responses:
                return self._responses[variable]
            strand = kwargs.get("strand")
            if strand in self._responses:
                return self._responses[strand]
            raise KeyError(f"No response configured for {kwargs}")
        return self._responses


class _RecordingQuantNado:
    def __init__(
        self,
        coverage: xr.DataArray,
        *,
        coverage_fwd: xr.DataArray | None = None,
        coverage_rev: xr.DataArray | None = None,
        methylation: xr.DataArray | None = None,
        variant_ref: xr.DataArray | None = None,
        variant_alt: xr.DataArray | None = None,
        variant_gt: xr.DataArray | None = None,
    ) -> None:
        self.coverage_calls: list[dict] = []
        self._coverage = coverage
        coverage_responses = {"+": coverage_fwd, "-": coverage_rev}
        self.coverage = _RecordingStore({k: v for k, v in coverage_responses.items() if v is not None})
        self.methylation = _RecordingStore(
            {"methylation_pct": methylation} if methylation is not None else {}
        )
        variant_responses = {}
        if variant_ref is not None:
            variant_responses["allele_depth_ref"] = variant_ref
        if variant_alt is not None:
            variant_responses["allele_depth_alt"] = variant_alt
        if variant_gt is not None:
            variant_responses["genotype"] = variant_gt
        self.variants = _RecordingStore(variant_responses)

    def extract_region(self, region: str, samples: list[str] | None = None, **kwargs):
        self.coverage_calls.append({"region": region, "samples": samples, **kwargs})
        return self._coverage


class TestQuantNadoTrackValidation:
    def test_coverage_requires_source(self):
        with pytest.raises(ValueError, match="requires one source"):
            QuantNadoCoverageTrack(sample="s1")

    def test_stranded_requires_source(self):
        with pytest.raises(ValueError, match="requires quantnado/dataset_path"):
            QuantNadoStrandedCoverageTrack(sample="s1")

    def test_methylation_requires_source(self):
        with pytest.raises(ValueError, match="requires one source"):
            QuantNadoMethylationTrack(sample="s1")

    def test_variant_requires_source(self):
        with pytest.raises(ValueError, match="requires quantnado/dataset_path"):
            QuantNadoVariantTrack(sample="s1")


class TestQuantNadoArrayMode:
    def test_coverage_plot_array_mode(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=105)
        track = QuantNadoCoverageTrack(
            sample="s1",
            coverage_data=_make_dataarray([0, 2, 4, 1, 3]),
            title="cov",
        )
        fig, ax = plt.subplots()
        track.plot(ax, gr)
        y_min, y_max = ax.get_ylim()
        assert y_min <= 0
        assert y_max >= 4
        assert len(ax.lines) >= 1
        plt.close(fig)

    def test_stranded_plot_array_mode(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=104)
        track = QuantNadoStrandedCoverageTrack(
            sample="s1",
            coverage_fwd_data=_make_dataarray([1, 3, 2, 4]),
            coverage_rev_data=_make_dataarray([2, 1, 3, 2]),
            title="strand",
        )
        fig, ax = plt.subplots()
        track.plot(ax, gr)
        y_min, y_max = ax.get_ylim()
        assert y_min < 0
        assert y_max > 0
        assert len(ax.lines) >= 2
        plt.close(fig)

    def test_methylation_plot_array_mode(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=103)
        track = QuantNadoMethylationTrack(
            sample="s1",
            methylation_data=_make_dataarray([10, 50, 90]),
            title="meth",
        )
        fig, ax = plt.subplots()
        track.plot(ax, gr)
        assert ax.get_ylim() == (0.0, 100.0)
        assert len(ax.collections) >= 1
        plt.close(fig)

    def test_variant_plot_array_mode(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=103)
        track = QuantNadoVariantTrack(
            sample="s1",
            allele_depth_ref_data=_make_dataarray([10, 10, 10]),
            allele_depth_alt_data=_make_dataarray([0, 5, 15]),
            title="var",
        )
        fig, ax = plt.subplots()
        track.plot(ax, gr)
        y_min, y_max = ax.get_ylim()
        assert y_min <= 0
        assert y_max >= 1
        assert len(ax.collections) >= 1
        plt.close(fig)


class TestQuantNadoObjectMode:
    def test_coverage_calls_quantnado_extract_region(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=105)
        da = _make_dataarray([1, 2, 3, 4, 5])
        qn = _RecordingQuantNado(coverage=da)
        track = QuantNadoCoverageTrack(sample="s1", quantnado=qn)
        fetched = track.fetch_data(gr)

        assert fetched["position"].tolist() == [100, 101, 102, 103, 104]
        assert qn.coverage_calls == [{"region": "chr1:100-105", "samples": ["s1"]}]

    def test_stranded_calls_coverage_store_with_strands(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=103)
        qn = _RecordingQuantNado(
            coverage=_make_dataarray([0, 0, 0]),
            coverage_fwd=_make_dataarray([1, 2, 3]),
            coverage_rev=_make_dataarray([3, 2, 1]),
        )
        track = QuantNadoStrandedCoverageTrack(sample="s1", quantnado=qn)
        fetched = track.fetch_data(gr)

        assert fetched["forward"].tolist() == [1.0, 2.0, 3.0]
        assert fetched["reverse"].tolist() == [3.0, 2.0, 1.0]
        assert qn.coverage.calls == [
            {"region": "chr1:100-103", "samples": ["s1"], "strand": "+"},
            {"region": "chr1:100-103", "samples": ["s1"], "strand": "-"},
        ]

    def test_methylation_calls_store_with_variable(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=103)
        qn = _RecordingQuantNado(
            coverage=_make_dataarray([0, 0, 0]),
            methylation=_make_dataarray([5, 15, 25]),
        )
        track = QuantNadoMethylationTrack(sample="s1", quantnado=qn)
        track.fetch_data(gr)

        assert qn.methylation.calls == [
            {
                "region": "chr1:100-103",
                "samples": ["s1"],
                "variable": "methylation_pct",
            }
        ]

    def test_variant_calls_store_with_expected_variables(self):
        gr = GenomicRegion(chromosome="chr1", start=100, end=103)
        qn = _RecordingQuantNado(
            coverage=_make_dataarray([0, 0, 0]),
            variant_ref=_make_dataarray([10, 10, 10]),
            variant_alt=_make_dataarray([2, 5, 8]),
            variant_gt=_make_dataarray([0, 1, 2]),
        )
        track = QuantNadoVariantTrack(sample="s1", quantnado=qn)
        track.fetch_data(gr)

        assert qn.variants.calls == [
            {
                "region": "chr1:100-103",
                "samples": ["s1"],
                "variable": "allele_depth_ref",
            },
            {
                "region": "chr1:100-103",
                "samples": ["s1"],
                "variable": "allele_depth_alt",
            },
            {
                "region": "chr1:100-103",
                "samples": ["s1"],
                "variable": "genotype",
            },
        ]

    def test_variant_fallback_genotype_classification(self):
        track = QuantNadoVariantTrack(
            sample="s1",
            allele_depth_ref_data=_make_dataarray([10, 10, 10]),
            allele_depth_alt_data=_make_dataarray([0, 6, 50]),
        )
        gr = GenomicRegion(chromosome="chr1", start=100, end=103)
        fetched = track.fetch_data(gr)
        total = fetched["ref_depth"] + fetched["alt_depth"]
        af = np.where(total > 0, fetched["alt_depth"] / total, np.nan)
        genotype = track._classify_genotype(af, None)

        assert genotype.tolist() == [-1, 1, 2]
