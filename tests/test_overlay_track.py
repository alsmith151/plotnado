import pandas as pd

from plotnado import GenomicFigure
from plotnado.tracks.base import GenomicRegion
from plotnado.tracks import BigWigTrack, BigwigOverlay, OverlayTrack
from plotnado.tracks.enums import PlotStyle


def _signal(scale: float = 1.0) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [1000, 1100, 1200],
            "end": [1100, 1200, 1300],
            "value": [1.0 * scale, 2.0 * scale, 1.5 * scale],
        }
    )


def test_overlay_track_initializes_private_instances():
    overlay = OverlayTrack(
        tracks=[
            BigWigTrack(data=_signal(1.0), title="a"),
            BigWigTrack(data=_signal(2.0), title="b"),
        ]
    )

    assert len(overlay._track_instances) == 2


def test_bigwigoverlay_backwards_compatible():
    overlay = BigwigOverlay(
        tracks=[
            BigWigTrack(data=_signal(1.0), title="a"),
            BigWigTrack(data=_signal(2.0), title="b"),
        ]
    )

    assert isinstance(overlay, OverlayTrack)


def test_figure_overlay_alias_plots():
    fig = GenomicFigure()
    fig.add_track(
        "overlay",
        tracks=[
            BigWigTrack(data=_signal(1.0), title="a"),
            BigWigTrack(data=_signal(2.0), title="b"),
        ],
        title="overlay",
    )

    out = fig.plot("chr1:1000-1300", show=False)
    assert out is not None


def test_overlay_shared_limits_include_all_components():
    overlay = OverlayTrack(
        tracks=[
            BigWigTrack(data=_signal(1.0), title="a"),
            BigWigTrack(data=_signal(10.0), title="b"),
        ]
    )

    limits = overlay._shared_limits(GenomicRegion(chromosome="chr1", start=1000, end=1300))

    assert limits == (0.0, 20.0)


def test_overlay_shared_limits_preserve_explicit_overrides():
    overlay = OverlayTrack(
        tracks=[
            BigWigTrack(data=_signal(1.0), title="a"),
            BigWigTrack(data=_signal(10.0), title="b"),
        ],
        aesthetics={"min_value": -5.0, "max_value": 7.0},
    )

    limits = overlay._shared_limits(GenomicRegion(chromosome="chr1", start=1000, end=1300))

    assert limits == (-5.0, 7.0)


def test_overlay_propagates_style_to_wrapped_bigwig_inputs():
    overlay = OverlayTrack(
        tracks=["track-a.bw", "track-b.bw"],
        aesthetics={"style": PlotStyle.FRAGMENT},
    )

    assert all(isinstance(track, BigWigTrack) for track in overlay._track_instances)
    assert all(track.style == PlotStyle.FRAGMENT for track in overlay._track_instances)
