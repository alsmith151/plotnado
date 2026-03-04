from __future__ import annotations

from pathlib import Path

import pandas as pd

from plotnado import Figure


def bed_intervals() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "start": [1_008_000, 1_020_000, 1_050_000, 1_066_000],
            "end": [1_014_000, 1_032_000, 1_061_000, 1_074_000],
            "name": ["a", "b", "c", "d"],
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


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = Figure(track_height=1.25)
    fig.scalebar()
    fig.axis()
    fig.bed(bed_intervals(), title="BED intervals", display="expanded", show_labels=True)
    fig.narrowpeak(
        narrowpeaks(),
        title="narrowPeak",
        color_by="signalValue",
        cmap="Oranges",
        show_labels=True,
        show_summit=True,
    )

    outfile = outdir / "track_bed_and_narrowpeak.png"
    fig.save(outfile, "chr1:1,005,000-1,078,000", dpi=220)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
