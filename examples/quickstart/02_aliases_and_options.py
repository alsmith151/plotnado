from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import Figure


def make_signal(phase: float = 0.0) -> pd.DataFrame:
    bins = np.arange(1_000_000, 1_080_000, 1_000)
    values = 4 + 1.5 * np.sin(np.linspace(phase, 5 + phase, bins.shape[0]))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + 1_000, "value": values})


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    print("Available aliases:")
    print(Figure.available_track_aliases())
    print("\nBigWig options (track/aesthetics/label):")
    print(Figure.track_options("bigwig"))

    fig = Figure().autocolor("Set2")
    fig.add_track("scalebar")
    fig.add_track("axis")
    fig.add_track("bigwig", data=make_signal(0.0), title="Replicate A", alpha=0.7, plot_title=True)
    fig.add_track("bigwig", data=make_signal(0.8), title="Replicate B", alpha=0.7, plot_title=True)

    outfile = outdir / "quickstart_aliases_and_options.png"
    fig.save(outfile, "chr1:1,005,000-1,070,000", dpi=200)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
