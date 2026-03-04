from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import Figure


def signal(style_shift: float) -> pd.DataFrame:
    bins = np.arange(1_000_000, 1_090_000, 1_000)
    values = 5 + 2.4 * np.sin(np.linspace(style_shift, 7 + style_shift, bins.shape[0]))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + 1_000, "value": values})


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = Figure(track_height=1.25)
    fig.scalebar()
    fig.bigwig(signal(0.0), title="fill", style="fill", color="#1f77b4")
    fig.bigwig(signal(0.8), title="fragment", style="fragment", color="#d62728")
    fig.bigwig(signal(1.6), title="scatter", style="scatter", color="#2ca02c", scatter_point_size=14)

    outfile = outdir / "track_bigwig_styles.png"
    fig.save(outfile, "chr1:1,010,000-1,080,000", dpi=220)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
