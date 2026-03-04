"""Basic PlotNado example with fully in-memory data."""

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import Figure


def synthetic_signal(start: int = 1_000_000, end: int = 1_100_000, step: int = 1_000) -> pd.DataFrame:
    bins = np.arange(start, end, step)
    values = 5.0 + 2.0 * np.sin(np.linspace(0, 6 * np.pi, bins.shape[0]))
    return pd.DataFrame(
        {
            "chrom": "chr1",
            "start": bins,
            "end": bins + step,
            "value": values,
        }
    )


def main() -> None:
    outdir = Path(__file__).parent / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = Figure(width=12, track_height=1.4)
    fig.scalebar(position="left")
    fig.axis()
    fig.genes("hg38", title="Genes")
    fig.bigwig(
        synthetic_signal(),
        title="Synthetic signal",
        style="fill",
        color="#1f77b4",
        alpha=0.8,
    )

    outfile = outdir / "basic_figure.png"
    fig.save(outfile, region="chr1:1,010,000-1,080,000", dpi=200)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
