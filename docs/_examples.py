from __future__ import annotations

import importlib.util

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plotnado import GenomicFigure

REGION = "chr1:1,010,000-1,080,000"
WIDE_REGION = "chr1:1,000,000-1,110,000"


def has_module(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def unavailable_figure(title: str, requirement: str):
    fig, ax = plt.subplots(figsize=(9, 1.8))
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        f"{title} requires {requirement}",
        ha="center",
        va="center",
        fontsize=10,
    )
    plt.close(fig)
    return fig


def signal(
    start: int = 1_000_000,
    end: int = 1_100_000,
    step: int = 1_000,
    phase: float = 0.0,
    scale: float = 1.0,
    baseline: float = 5.0,
) -> pd.DataFrame:
    bins = np.arange(start, end, step)
    values = scale * (
        baseline + 2.0 * np.sin(np.linspace(phase, 6 * np.pi + phase, bins.shape[0]))
    )
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + step, "value": values})


def review_signal(scale: float = 1.0, phase: float = 0.0) -> pd.DataFrame:
    bins = np.arange(1_000_000, 1_120_000, 1_000)
    values = scale * (1.2 + np.sin(np.linspace(phase, 6 + phase, bins.shape[0])))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + 1_000, "value": values})


def intervals() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "start": [1_008_000, 1_020_000, 1_050_000, 1_066_000],
            "end": [1_014_000, 1_032_000, 1_061_000, 1_074_000],
            "name": ["enhancer_a", "enhancer_b", "promoter", "domain"],
        }
    )


def narrowpeaks() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [1_012_000, 1_038_000, 1_060_000],
            "end": [1_018_000, 1_047_000, 1_070_000],
            "name": ["np1", "np2", "np3"],
            "score": [300, 700, 500],
            "strand": [".", ".", "."],
            "signalValue": [12.0, 48.0, 30.0],
            "pValue": [5.2, 12.3, 8.1],
            "qValue": [4.1, 10.0, 6.2],
            "peak": [1200, 1800, 2200],
        }
    )


def links() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom1": ["chr1", "chr1", "chr1"],
            "start1": [1_010_000, 1_022_000, 1_042_000],
            "end1": [1_012_000, 1_024_000, 1_045_000],
            "chrom2": ["chr1", "chr1", "chr1"],
            "start2": [1_035_000, 1_054_000, 1_072_000],
            "end2": [1_037_000, 1_056_000, 1_074_000],
            "score": [2.2, 6.5, 9.8],
        }
    )


def quickstart_figure() -> GenomicFigure:
    fig = GenomicFigure(width=11, track_height=1.25)
    fig.scalebar()
    fig.axis()
    fig.bigwig(signal(scale=1.15), title="Synthetic signal", style="fill", color="#1f77b4")
    fig.bed(intervals(), title="Intervals", display="expanded", show_labels=True)
    return fig


def style_comparison() -> GenomicFigure:
    fig = GenomicFigure(track_height=1.15)
    fig.scalebar()
    fig.bigwig(signal(phase=0.0), title="fill", style="fill", color="#1f77b4")
    fig.bigwig(signal(phase=0.8), title="fragment", style="fragment", color="#d62728")
    fig.bigwig(
        signal(phase=1.6),
        title="scatter",
        style="scatter",
        color="#2ca02c",
        scatter_point_size=10,
    )
    fig.bigwig(signal(phase=2.4), title="std", style="std", color="#9467bd")
    return fig


def overlay_comparison() -> GenomicFigure:
    fig = GenomicFigure(track_height=1.2)
    fig.autoscale(True)
    fig.highlight("chr1:1,032,000-1,046,000")
    fig.highlight_style(color="#ffdd57", alpha=0.22)
    fig.axis()
    fig.bigwig(review_signal(2.0), title="Control", autoscale_group="signal", color="#1f77b4")
    fig.bigwig(
        review_signal(10.0, 1.2),
        title="Treatment",
        autoscale_group="signal",
        color="#d62728",
    )
    fig.overlay(
        [review_signal(5.5, 2.0), review_signal(6.5, 2.8)],
        title="Overlay",
        autoscale_group="signal",
        colors=["#2ca02c", "#9467bd"],
        alpha=0.55,
    )
    return fig
