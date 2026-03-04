# Data Inputs

## In-memory tables

Many tracks accept `pandas.DataFrame` inputs directly.

### BigWig-like signal

Expected columns:

- `chrom`
- `start`
- `end`
- `value`

### BED-like intervals

Expected columns (minimum):

- `chrom`
- `start`
- `end`

Optional labels use `name` (or a custom `label_field`).

### Links/BEDPE-like interactions

Expected columns:

- `chrom1`, `start1`, `end1`
- `chrom2`, `start2`, `end2`
- optional `score`

## File-based inputs

- BigWig tracks: BigWig files (`.bw`/`.bigWig`)
- BED tracks: BED/BigBed
- narrowPeak tracks: `.narrowPeak`
- Matrix tracks: Cooler/MCool

Use in-memory data for quick prototyping and tests, then switch to file-backed inputs for real datasets.
