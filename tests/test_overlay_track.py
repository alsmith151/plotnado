import pandas as pd

from plotnado import GenomicFigure
from plotnado.tracks import BigWigTrack, BigwigOverlay, OverlayTrack


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
