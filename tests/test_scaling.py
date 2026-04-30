import matplotlib.axes
import numpy as np
import pandas as pd
from pydantic import BaseModel, Field

from plotnado.tracks import BigWigTrack, OverlayTrack
from plotnado.tracks.base import GenomicRegion, Track
from plotnado.tracks.scaling import Autoscaler, calculate_data_limits


class _AutoscaleAesthetics(BaseModel):
    min_value: float | None = None
    max_value: float | None = None


class _AutoscaleTrack(Track):
    aesthetics: _AutoscaleAesthetics = Field(default_factory=_AutoscaleAesthetics)

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:  # noqa: ARG002
        return pd.DataFrame({"value": [1.0, 2.5, 4.0]})

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:  # noqa: ARG002
        return None


class _NonNumericTrack(Track):
    aesthetics: _AutoscaleAesthetics = Field(default_factory=_AutoscaleAesthetics)

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:  # noqa: ARG002
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "meta": [[1, 2], [3, 4]],
            }
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:  # noqa: ARG002
        return None


class _NoAutoscaleAestheticTrack(Track):
    class Aesthetics(BaseModel):
        color: str = "black"

    aesthetics: Aesthetics = Field(default_factory=Aesthetics)

    def fetch_data(self, gr: GenomicRegion) -> np.ndarray:  # noqa: ARG002
        return np.array([1000.0, 2000.0, 3000.0])

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:  # noqa: ARG002
        return None


def test_autoscaler_ignores_non_numeric_dataframe_columns() -> None:
    gr = GenomicRegion(chromosome="chr1", start=100, end=300)
    numeric_track = _AutoscaleTrack()
    non_numeric_track = _NonNumericTrack()

    scaler = Autoscaler(tracks=[numeric_track, non_numeric_track], gr=gr)
    scaler.apply()

    assert numeric_track.aesthetics.min_value == 0.0
    assert numeric_track.aesthetics.max_value == 4.0
    assert non_numeric_track.aesthetics.min_value == 0.0
    assert non_numeric_track.aesthetics.max_value == 4.0


def test_autoscaler_uses_only_tracks_with_min_max_aesthetics() -> None:
    gr = GenomicRegion(chromosome="chr1", start=100, end=300)
    numeric_track = _AutoscaleTrack()
    outlier_track = _NoAutoscaleAestheticTrack()

    scaler = Autoscaler(tracks=[numeric_track, outlier_track], gr=gr)
    scaler.apply()

    assert numeric_track.aesthetics.min_value == 0.0
    assert numeric_track.aesthetics.max_value == 4.0


def test_calculate_data_limits_flattens_overlay_like_payloads() -> None:
    overlay_payload = [
        pd.DataFrame({"value": [1.0, 2.5, 4.0]}),
        pd.DataFrame({"value": [10.0, 12.0]}),
    ]

    assert calculate_data_limits(overlay_payload) == (0.0, 12.0)


def test_autoscaler_aggregates_overlay_component_values() -> None:
    gr = GenomicRegion(chromosome="chr1", start=100, end=300)
    base_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [200, 300],
            "value": [2.0, 4.0],
        }
    )
    overlay_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [200, 300],
            "value": [10.0, 20.0],
        }
    )

    numeric_track = BigWigTrack(data=base_df)
    overlay_track = OverlayTrack(
        tracks=[
            BigWigTrack(data=overlay_df, title="overlay-a"),
            BigWigTrack(data=overlay_df * 0 + overlay_df.assign(value=[5.0, 8.0]), title="overlay-b"),
        ]
    )

    scaler = Autoscaler(tracks=[numeric_track, overlay_track], gr=gr)
    scaler.apply()

    assert numeric_track.aesthetics.min_value == 0.0
    assert numeric_track.aesthetics.max_value == 20.0
    assert overlay_track.aesthetics.min_value == 0.0
    assert overlay_track.aesthetics.max_value == 20.0


def test_autoscaler_preserves_explicit_overlay_limits() -> None:
    gr = GenomicRegion(chromosome="chr1", start=100, end=300)
    regular_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [200, 300],
            "value": [2.0, 4.0],
        }
    )
    overlay_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [200, 300],
            "value": [10.0, 20.0],
        }
    )

    numeric_track = BigWigTrack(data=regular_df)
    overlay_track = OverlayTrack(
        tracks=[BigWigTrack(data=overlay_df, title="overlay")],
        aesthetics={"min_value": -5.0, "max_value": 7.0},
    )

    scaler = Autoscaler(tracks=[numeric_track, overlay_track], gr=gr)
    scaler.apply()

    assert numeric_track.aesthetics.min_value == 0.0
    assert numeric_track.aesthetics.max_value == 20.0
    assert overlay_track.aesthetics.min_value == -5.0
    assert overlay_track.aesthetics.max_value == 7.0
