from __future__ import annotations

from pathlib import Path
import random
import pandas as pd
from plotnado import Figure


def dense_gene_annotations() -> pd.DataFrame:
    records = []
    starts = [1_010_000, 1_013_500, 1_017_000, 1_020_500, 1_024_000, 1_028_000, 1_031_500]
    lengths = [6_000, 5_500, 6_300, 5_200, 5_800, 5_400, 6_100]

    for index, (start, length) in enumerate(zip(starts, lengths), start=1):
        end = start + length
        records.append(
            {
                "chrom": "chr1",
                "start": start,
                "end": end,
                "geneid": f"GENE_{index}",
                "strand": "+" if index % 2 else "-",
                "block_count": 3,
                "block_starts": [0, int(length * 0.42), int(length * 0.78)],
                "block_sizes": [int(length * 0.16), int(length * 0.12), int(length * 0.18)],
            }
        )

    return pd.DataFrame.from_records(records)


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    genes_df = dense_gene_annotations()

    fig = Figure(width=12, track_height=1.2)
    fig.scalebar()
    fig.axis()

    fig.genes(
        data=genes_df,
        title="Method: Staggered labels + connectors",
        display="collapsed",
        label_overlap_strategy="stagger",
        label_max_chars=120,
    )
    fig.genes(
        data=genes_df,
        title="Method: Suppress overlaps (offset labels)",
        display="collapsed",
        label_overlap_strategy="suppress",
        label_max_chars=120,
    )
    fig.genes(
        data=genes_df,
        title="Method: Auto-expand rows + suppress overlaps",
        display="collapsed",
        label_overlap_strategy="auto_expand",
        max_number_of_rows=6,
        label_max_chars=120,
    )

    region = "chr1:1,006,000-1,040,000"
    outfile = outdir / "recipe_gene_label_strategies.png"
    fig.save(outfile, region=region, dpi=220)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()
