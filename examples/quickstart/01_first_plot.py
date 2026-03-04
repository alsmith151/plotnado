from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import Figure


def synthetic_signal(start: int = 1_000_000, end: int = 1_100_000, step: int = 1_000) -> pd.DataFrame:
    bins = np.arange(start, end, step)
    values = 6 + 2 * np.sin(np.linspace(0, 6 * np.pi, bins.shape[0]))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + step, "value": values})


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = Figure(width=11)
    fig.scalebar()
    fig.axis()
    fig.genes("hg38", title="Genes")
    fig.bigwig(synthetic_signal(), title="Synthetic ChIP", color="#1f77b4", style="fill")

    outfile = outdir / "quickstart_first_plot.png"
    fig.save(outfile, "chr1:1,010,000-1,080,000", dpi=200)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
