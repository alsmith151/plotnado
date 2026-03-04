from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotnado import Figure, Theme


def signal() -> pd.DataFrame:
    bins = np.arange(1_000_000, 1_090_000, 1_000)
    values = 4 + 2 * np.cos(np.linspace(0, 5, bins.shape[0]))
    return pd.DataFrame({"chrom": "chr1", "start": bins, "end": bins + 1_000, "value": values})


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = Figure(theme=Theme.publication())
    fig.scalebar(plot_scale=True, scale_size=10)
    fig.axis(num_ticks=4)
    fig.bigwig(
        signal(),
        title="Publication theme",
        color="#2c7fb8",
        label_box_enabled=True,
        label_box_alpha=0.85,
        title_location="left",
        scale_location="right",
    )

    image_out = outdir / "recipe_theme_labels.png"
    fig.save(image_out, "chr1:1,010,000-1,080,000", dpi=220)
    print(f"Saved {image_out}")

    toml_fig = Figure(theme=Theme.publication())
    toml_fig.scalebar(plot_scale=True, scale_size=10)
    toml_fig.axis(num_ticks=4)
    toml_fig.genes("hg38", title="Genes", title_location="left")

    toml_out = outdir / "recipe_theme_labels.toml"
    try:
        toml_fig.to_toml(str(toml_out))
        loaded = Figure.from_toml(str(toml_out))
        reloaded_out = outdir / "recipe_theme_labels_from_toml.png"
        loaded.save(reloaded_out, "chr1:1,010,000-1,080,000", dpi=220)
        print(f"Saved {toml_out} and {reloaded_out}")
    except ImportError:
        print("Skipping TOML roundtrip. Install optional dependency 'tomli-w' to enable it.")


if __name__ == "__main__":
    main()
