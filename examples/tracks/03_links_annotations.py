from __future__ import annotations

from pathlib import Path

import pandas as pd

from plotnado import GenomicFigure


def link_df() -> pd.DataFrame:
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


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    fig = GenomicFigure(track_height=1.2)
    fig.axis()
    fig.links(link_df(), title="Loops", color_by_score=True, cmap="viridis", alpha=0.8)
    fig.hline(0.1, color="#666666", linewidth=0.8)
    fig.vline(1_048_000, color="#222222", linewidth=1.0)

    outfile = outdir / "track_links_annotations.png"
    fig.save(outfile, "chr1:1,005,000-1,078,000", dpi=220)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
