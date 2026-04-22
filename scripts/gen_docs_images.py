"""Generate all images used in the PlotNado documentation.

Run with:  uv run python scripts/gen_docs_images.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import GenomicFigure

REGION = "chr1:1,010,000-1,080,000"
DPI = 180
OUT = Path(__file__).resolve().parents[1] / "docs" / "images"


def _signal(phase: float = 0.0, scale: float = 1.0, step: int = 1_000) -> pd.DataFrame:
    bins = np.arange(1_000_000, 1_100_000, step)
    values = scale * (4 + 2 * np.sin(np.linspace(phase, 6 + phase, bins.shape[0])))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + step, "value": values})


def _bed() -> pd.DataFrame:
    return pd.DataFrame({
        "chrom": ["chr1"] * 5,
        "start": [1_012_000, 1_022_000, 1_038_000, 1_055_000, 1_068_000],
        "end":   [1_016_000, 1_028_000, 1_046_000, 1_062_000, 1_074_000],
        "name":  ["a", "b", "c", "d", "e"],
    })


def _narrowpeaks() -> pd.DataFrame:
    return pd.DataFrame({
        "chrom":       ["chr1", "chr1", "chr1"],
        "start":       [1_015_000, 1_040_000, 1_060_000],
        "end":         [1_022_000, 1_050_000, 1_070_000],
        "name":        ["np1", "np2", "np3"],
        "score":       [300, 700, 500],
        "strand":      [".", ".", "."],
        "signalValue": [12.0, 48.0, 30.0],
        "pValue":      [5.2, 12.3, 8.1],
        "qValue":      [4.1, 10.0, 6.2],
        "peak":        [1200, 1800, 2200],
    })


# ── helpers ──────────────────────────────────────────────────────────────────

def save(fig: GenomicFigure, name: str, subdir: str = "") -> None:
    folder = OUT / subdir if subdir else OUT
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    fig.save(path, REGION, dpi=DPI)
    print(f"  {path.relative_to(OUT.parent.parent)}")


# ── quickstart image ──────────────────────────────────────────────────────────

def gen_quickstart() -> None:
    print("quickstart...")
    fig = GenomicFigure(width=11)
    fig.scalebar()
    fig.axis()
    fig.genes("hg38", title="Genes")
    fig.bigwig(_signal(), title="ChIP signal", color="#1f77b4", style="fill")
    save(fig, "quickstart.png")


# ── bigwig styles ─────────────────────────────────────────────────────────────

def gen_bigwig_styles() -> None:
    print("bigwig styles...")
    fig = GenomicFigure(track_height=1.25)
    fig.scalebar()
    fig.bigwig(_signal(0.0), title="fill",     style="fill",     color="#1f77b4")
    fig.bigwig(_signal(0.8), title="fragment", style="fragment", color="#d62728")
    fig.bigwig(_signal(1.6), title="scatter",  style="scatter",  color="#2ca02c", scatter_point_size=10)
    fig.bigwig(_signal(2.4), title="std",      style="std",      color="#9467bd")
    save(fig, "bigwig_styles.png", "aesthetics")


# ── color and alpha ───────────────────────────────────────────────────────────

def gen_color_alpha() -> None:
    print("color and alpha...")
    fig = GenomicFigure(track_height=1.1)
    fig.scalebar()
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e"]
    alphas = [1.0, 0.8, 0.5, 0.2]
    for i, (color, alpha) in enumerate(zip(colors, alphas)):
        fig.bigwig(_signal(i * 0.5), title=f"color={color}  alpha={alpha}", style="fill", color=color, alpha=alpha)
    save(fig, "color_alpha.png", "aesthetics")


# ── label placement ───────────────────────────────────────────────────────────

def gen_label_placement() -> None:
    print("label placement...")
    fig = GenomicFigure(track_height=1.2)
    fig.scalebar()
    configs = [
        dict(title_location="left",   label_on_track=False, title="title_location='left'"),
        dict(title_location="right",  label_on_track=False, title="title_location='right'"),
        dict(title_location="center", label_on_track=False, title="title_location='center'"),
        dict(title_location="left",   label_on_track=True,  title="label_on_track=True"),
    ]
    for i, cfg in enumerate(configs):
        fig.bigwig(_signal(i * 0.6), style="fill", color="#1f77b4", **cfg)
    save(fig, "label_placement.png", "aesthetics")


# ── scale / label box ─────────────────────────────────────────────────────────

def gen_scale_box() -> None:
    print("scale and label box...")
    fig = GenomicFigure(track_height=1.3)
    fig.scalebar(plot_scale=True, scale_size=10)
    fig.bigwig(
        _signal(),
        title="label_box + scale (right)",
        color="#2c7fb8",
        label_box_enabled=True,
        label_box_alpha=0.85,
        title_location="left",
        plot_scale=True,
        scale_location="right",
    )
    fig.bigwig(
        _signal(1.5),
        title="no box, scale left",
        color="#d62728",
        label_box_enabled=False,
        plot_scale=True,
        scale_location="left",
    )
    save(fig, "scale_box.png", "aesthetics")


# ── BED & narrowPeak ──────────────────────────────────────────────────────────

def gen_bed_narrowpeak() -> None:
    print("BED & narrowPeak...")
    fig = GenomicFigure(track_height=1.25)
    fig.scalebar()
    fig.axis()
    fig.bed(_bed(), title="BED  (show_labels=True)", display="expanded", show_labels=True, color="#4dac26")
    fig.narrowpeak(
        _narrowpeaks(),
        title="narrowPeak  (color_by='signalValue')",
        color_by="signalValue",
        cmap="Oranges",
        show_labels=True,
        show_summit=True,
    )
    save(fig, "bed_narrowpeak.png", "aesthetics")


# ── overlay + autoscale + highlight ──────────────────────────────────────────

def gen_overlay() -> None:
    print("overlay / autoscale / highlight...")
    fig = GenomicFigure(track_height=1.25)
    fig.autoscale(True)
    fig.highlight("chr1:1,032,000-1,046,000")
    fig.highlight_style(color="#ffdd57", alpha=0.22)
    fig.axis()
    fig.bigwig(_signal(0.0, scale=2.0),  title="Sample A (autoscale_group='g1')", autoscale_group="g1", style="fill", color="#1f77b4")
    fig.bigwig(_signal(1.2, scale=10.0), title="Sample B (autoscale_group='g1')", autoscale_group="g1", style="fill", color="#d62728")
    fig.overlay(
        [_signal(0.9, scale=4.2), _signal(1.5, scale=3.5)],
        title="Overlay",
        colors=["#2ca02c", "#9467bd"],
        alpha=0.55,
    )
    save(fig, "overlay_autoscale.png", "aesthetics")


# ── autocolor + color_group ───────────────────────────────────────────────────

def gen_autocolor() -> None:
    print("autocolor / color_group...")
    fig = GenomicFigure(track_height=1.1)
    fig.autocolor()
    fig.scalebar()
    fig.bigwig(_signal(0.0), title="Sample A signal", color_group="A", style="fill")
    fig.bed(_bed(), title="Sample A peaks", color_group="A", display="expanded")
    fig.bigwig(_signal(1.5), title="Sample B signal", color_group="B", style="fill")
    fig.bed(_bed(), title="Sample B peaks", color_group="B", display="expanded")
    save(fig, "autocolor.png", "aesthetics")


# ── gene track display modes ──────────────────────────────────────────────────

def gen_gene_display() -> None:
    print("gene display modes...")
    fig = GenomicFigure(track_height=1.6)
    fig.scalebar()
    fig.genes("hg38", title="display='collapsed'", display="collapsed")
    fig.genes("hg38", title="display='expanded'",  display="expanded")
    save(fig, "gene_display.png", "aesthetics")


# ── links ─────────────────────────────────────────────────────────────────────

def gen_links() -> None:
    print("links / hline / vline...")
    links_df = pd.DataFrame({
        "chrom1":  ["chr1", "chr1", "chr1"],
        "start1":  [1_010_000, 1_022_000, 1_042_000],
        "end1":    [1_012_000, 1_024_000, 1_045_000],
        "chrom2":  ["chr1", "chr1", "chr1"],
        "start2":  [1_035_000, 1_054_000, 1_072_000],
        "end2":    [1_037_000, 1_056_000, 1_074_000],
        "score":   [2.2, 6.5, 9.8],
    })
    fig = GenomicFigure(track_height=1.2)
    fig.axis()
    fig.links(links_df, title="Links (color_by_score=True)", color_by_score=True, cmap="viridis", alpha=0.8)
    fig.hline(0.1, color="#666666", linewidth=0.8)
    fig.vline(1_048_000, color="#222222", linewidth=1.0)
    save(fig, "links.png", "aesthetics")


if __name__ == "__main__":
    print("Generating documentation images...")
    gen_quickstart()
    gen_bigwig_styles()
    gen_color_alpha()
    gen_label_placement()
    gen_scale_box()
    gen_bed_narrowpeak()
    gen_overlay()
    gen_autocolor()
    gen_gene_display()
    gen_links()
    print("Done.")
