"""Advanced PlotNado features with runnable synthetic data."""

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import Figure


def make_signal(scale: float, phase: float = 0.0) -> pd.DataFrame:
    bins = np.arange(1_000_000, 1_120_000, 1_000)
    signal = scale * (1.2 + np.sin(np.linspace(phase, 6 + phase, bins.shape[0])))
    return pd.DataFrame(
        {
            "chrom": "chr1",
            "start": bins,
            "end": bins + 1_000,
            "value": signal,
        }
    )


def make_intervals() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "start": [1_008_000, 1_025_000, 1_055_000, 1_085_000],
            "end": [1_017_000, 1_040_000, 1_067_000, 1_098_000],
            "name": ["peak_a", "peak_b", "peak_c", "peak_d"],
        }
    )


def main() -> None:
    outdir = Path(__file__).parent / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = Figure(width=12, track_height=1.3)
    fig.autocolor("tab10")
    fig.highlight("chr1:1,030,000-1,045,000")
    fig.highlight_style(color="#f6d55c", alpha=0.2)

    fig.scalebar(position="right")
    fig.axis(num_ticks=5)
    fig.genes("hg38", title="Genes")
    fig.bed(make_intervals(), title="Intervals", display="expanded", show_labels=True)
    fig.bigwig(make_signal(4.0), title="Control", style="fill", autoscale_group="signal")
    fig.bigwig(make_signal(5.2, phase=0.8), title="Treatment", style="fill", autoscale_group="signal")
    fig.bigwig(make_signal(3.0, phase=1.3), title="Input", style="fragment")
    fig.vline(1_060_000, color="#333333", linewidth=1.0)
    fig.hline(0.0, color="#666666", linewidth=0.6)

    outfile = outdir / "advanced_features.png"
    fig.save(outfile, region="chr1:1,000,000-1,110,000", dpi=220)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
