"""
Synthetic and example data helpers for plotnado documentation and tutorials.

These functions return small in-memory DataFrames that mimic real genomic
data formats. Swap them for real files or URLs in production:

    signal()      → BigWig path / URL, or DataFrame(chrom, start, end, value)
    intervals()   → BED path / URL, or DataFrame(chrom, start, end, name)
    narrowpeaks() → .narrowPeak path, or DataFrame with ENCODE narrowPeak columns
    links()       → BEDPE path, or DataFrame(chrom1, start1, end1, chrom2, start2, end2, score)
"""

from __future__ import annotations

import importlib.util
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import GenomicFigure, PlotStyle

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REGION = "chr1:1,010,000-1,080,000"
WIDE_REGION = "chr1:1,000,000-1,110,000"

BLUEPRINT_MONOCYTE_BW = (
    "http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood"
    "/C000S5/CD14-positive_CD16-negative_classical_monocyte/ChIP-Seq/NCMLS"
    "/C000S5H2.ERX173536.H3K27ac.bwa.GRCh38.20150529.bw"
)
BLUEPRINT_MONOCYTE_BB = (
    "http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood"
    "/C000S5/CD14-positive_CD16-negative_classical_monocyte/ChIP-Seq/NCMLS"
    "/C000S5H2.ERX173536.H3K27ac.bwa.GRCh38.20150527.bb"
)
BLUEPRINT_REGION = "chr12:69,310,000-69,400,000"

_CANDIDATE_BED = str(Path(__file__).parent / "data" / "blueprint_monocyte_lyz_candidate_peaks.bed")

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def has_module(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


# ---------------------------------------------------------------------------
# Synthetic signal generators
# ---------------------------------------------------------------------------


def _chip_signal(
    bins: np.ndarray,
    peak_centers: Sequence[int],
    peak_heights: Sequence[float],
    peak_widths: Sequence[float],
    noise_seed: int = 0,
    noise_smooth: int = 7,
) -> np.ndarray:
    """Sparse sharp Gaussian peaks on a smoothed exponential-noise background."""
    rng = np.random.default_rng(noise_seed)
    n = len(bins)
    noise = rng.exponential(0.25, n) + rng.exponential(0.05, n)
    # Smooth noise so adjacent bins are correlated (mimics genomic read pileup)
    if noise_smooth > 1:
        kernel = np.ones(noise_smooth) / noise_smooth
        noise = np.convolve(noise, kernel, mode="same")
    values = np.clip(noise, 0, None)
    for center, height, width in zip(peak_centers, peak_heights, peak_widths):
        values += height * np.exp(-0.5 * ((bins - center) / width) ** 2)
    return values


def signal(
    start: int = 1_000_000,
    end: int = 1_100_000,
    step: int = 200,
    phase: float = 0.0,
    scale: float = 1.0,
    baseline: float = 5.0,
) -> pd.DataFrame:
    """ChIP-seq-like synthetic signal over chr1.

    Returns DataFrame(chrom, start, end, value).
    Replace with a BigWig path/URL or your own DataFrame.
    """
    bins = np.arange(start, end, step)
    span = end - start
    peak_centers = [
        int(start + span * (0.20 + phase * 0.04)),
        int(start + span * (0.48 + phase * 0.03)),
        int(start + span * (0.73 + phase * 0.02)),
    ]
    peak_heights = [h * scale * baseline / 7 for h in [8.0, 14.0, 6.0]]
    peak_widths = [2_500, 3_000, 2_000]
    values = _chip_signal(bins, peak_centers, peak_heights, peak_widths, noise_seed=int(phase * 100))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + step, "value": values})


def review_signal(scale: float = 1.0, phase: float = 0.0) -> pd.DataFrame:
    """Wider ChIP-seq-like signal for overlay/autoscale examples.

    Returns DataFrame(chrom, start, end, value).
    Replace with a BigWig path/URL or your own DataFrame.
    """
    bins = np.arange(1_000_000, 1_120_000, 200)
    span = 120_000
    peak_centers = [
        int(1_000_000 + span * (0.22 + phase * 0.03)),
        int(1_000_000 + span * (0.55 + phase * 0.02)),
        int(1_000_000 + span * (0.78 + phase * 0.04)),
    ]
    peak_heights = [h * scale for h in [12.0, 20.0, 8.0]]
    peak_widths = [3_000, 4_000, 2_500]
    values = _chip_signal(bins, peak_centers, peak_heights, peak_widths, noise_seed=int(scale * 7 + phase * 13))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + 200, "value": values})


# ---------------------------------------------------------------------------
# Synthetic interval data
# ---------------------------------------------------------------------------


def intervals() -> pd.DataFrame:
    """Four synthetic genomic intervals on chr1.

    Returns DataFrame(chrom, start, end, name).
    Replace with a BED/BigBed path, URL, or your own DataFrame.
    """
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "start": [1_008_000, 1_020_000, 1_050_000, 1_066_000],
            "end": [1_014_000, 1_032_000, 1_061_000, 1_074_000],
            "name": ["enhancer_a", "enhancer_b", "promoter", "domain"],
        }
    )


def narrowpeaks() -> pd.DataFrame:
    """Three synthetic ENCODE narrowPeak rows on chr1.

    Returns DataFrame with columns: chrom, start, end, name, score, strand,
    signalValue, pValue, qValue, peak.
    Replace with a .narrowPeak path or your own DataFrame.
    """
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
    """Three synthetic interaction links on chr1.

    Returns DataFrame(chrom1, start1, end1, chrom2, start2, end2, score).
    Replace with a BEDPE path or your own DataFrame.
    """
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


# ---------------------------------------------------------------------------
# Pre-built figure helpers (used in docs; inline the code for notebooks)
# ---------------------------------------------------------------------------


def hero_figure() -> GenomicFigure:
    """Real Blueprint monocyte H3K27ac figure used as the homepage hero plot."""
    fig = GenomicFigure(width=11, track_height=0.7, theme="publication")
    fig.scalebar()
    fig.genes("hg38", height=0.55)
    fig.bigwig(
        BLUEPRINT_MONOCYTE_BW,
        title="Monocyte H3K27ac",
        style=PlotStyle.FRAGMENT,
        color="#d9485f",
        label_on_track=True,
        label_box_enabled=True,
        label_box_alpha=0.92,
        title_height=0.5,
        scale_height=0.5,
        plot_scale=True,
        height=0.75,
    )
    fig.bed(
        BLUEPRINT_MONOCYTE_BB,
        title="H3K27ac peaks",
        color="#f59e0b",
        draw_edges=False,
        height=0.35,
        label_on_track=True,
        label_box_enabled=True,
        label_box_alpha=0.92,
        title_height=0.7,
        show_labels=False,
    )
    fig.bed(
        _CANDIDATE_BED,
        title="Candidate loci",
        color="#0f766e",
        draw_edges=True,
        show_labels=True,
        label_field="name",
        font_size=7,
        height=0.4,
        label_on_track=True,
        label_box_enabled=True,
        label_box_alpha=0.92,
        title_height=0.7,
    )
    fig.axis()
    return fig


def quickstart_figure() -> GenomicFigure:
    fig = GenomicFigure(width=11, track_height=1.25)
    fig.scalebar()
    fig.bigwig(signal(scale=1.15), title="Synthetic signal", style="fill", color="#1f77b4")
    fig.bed(intervals(), title="Intervals", display="expanded", show_labels=True)
    fig.axis()
    return fig


def style_comparison() -> GenomicFigure:
    fig = GenomicFigure(track_height=1.15)
    fig.scalebar()
    fig.bigwig(signal(phase=0.0), title="fill", style="fill", color="#1f77b4")
    fig.bigwig(signal(phase=0.8), title="fragment", style="fragment", color="#d62728")
    fig.bigwig(signal(phase=1.6), title="scatter", style="scatter", color="#2ca02c", scatter_point_size=10)
    fig.bigwig(signal(phase=2.4), title="std", style="std", color="#9467bd")
    return fig


def overlay_comparison() -> GenomicFigure:
    fig = GenomicFigure(track_height=1.2)
    fig.autoscale(True)
    fig.highlight("chr1:1,032,000-1,046,000")
    fig.highlight_style(color="#ffdd57", alpha=0.22)
    fig.bigwig(review_signal(2.0), title="Control", autoscale_group="signal", color="#1f77b4")
    fig.bigwig(review_signal(10.0, 1.2), title="Treatment", autoscale_group="signal", color="#d62728")
    fig.overlay(
        [review_signal(5.5, 2.0), review_signal(6.5, 2.8)],
        title="Overlay",
        autoscale_group="signal",
        colors=["#2ca02c", "#9467bd"],
        alpha=0.55,
    )
    return fig
