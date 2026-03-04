from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import Figure


def signal(scale: float, phase: float = 0.0) -> pd.DataFrame:
    bins = np.arange(1_000_000, 1_110_000, 1_000)
    values = scale * (1.5 + np.sin(np.linspace(phase, 6 + phase, bins.shape[0])))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + 1_000, "value": values})


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = Figure(track_height=1.25)
    fig.autoscale(True)
    fig.highlight("chr1:1,032,000-1,046,000")
    fig.highlight_style(color="#ffdd57", alpha=0.22)

    fig.axis()
    fig.bigwig(signal(2.0), title="Group A", autoscale_group="g1", style="fill", color="#1f77b4")
    fig.bigwig(signal(10.0, 1.2), title="Group B", autoscale_group="g1", style="fill", color="#d62728")
    fig.overlay(
        [
            signal(4.2, 0.9),
            signal(3.5, 1.5),
        ],
        title="Overlay",
        colors=["#2ca02c", "#9467bd"],
        alpha=0.55,
    )

    outfile = outdir / "recipe_autoscale_overlay_highlight.png"
    fig.save(outfile, "chr1:1,005,000-1,085,000", dpi=220)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
