# Track Catalog

This page summarizes major track types and points to runnable scripts.

## Structural tracks

- `scalebar` (`ScaleBar`): genomic distance bar.
- `axis` (`GenomicAxis`): coordinate ticks and labels.
- `spacer` (`Spacer`): layout spacing.

## Signal tracks

- `bigwig` (`BigWigTrack`): continuous or interval signal (`style`: `fill`, `fragment`, `scatter`).
- `overlay` (`OverlayTrack`): generic overlay panel for multiple tracks/signals.
- `bigwig_overlay` (`BigwigOverlay`): backward-compatible alias of `OverlayTrack`.
- `bigwig_collection` (`BigWigCollection`): collection workflows for multiple BigWigs.
- `bigwig_diff` (`BigWigDiff`): subtraction/ratio/log2ratio between two tracks.

Runnable scripts:

- `examples/tracks/01_bigwig_styles.py`
- `examples/recipes/01_autoscale_overlay_highlight.py`

![BigWig styles](images/examples/track_bigwig_styles.png)
![Overlay recipe](images/examples/recipe_autoscale_overlay_highlight.png)

## Interval and annotation tracks

- `bed` (`BedTrack`): interval rectangles with optional labels.
- `narrowpeak` (`NarrowPeakTrack`): peak-specific intervals with summit and score-based colors.
- `genes` (`Genes`): gene models from bundled genome annotations or custom files.
- `links` (`LinksTrack`): arcs for BEDPE-style interactions.
- `hline` / `vline`: reference lines.
- `highlight` (`HighlightsFromFile`): region overlays.

Runnable scripts:

- `examples/tracks/02_bed_and_narrowpeak.py`
- `examples/tracks/03_links_annotations.py`

![BED and narrowPeak](images/examples/track_bed_and_narrowpeak.png)
![Links and annotations](images/examples/track_links_annotations.png)

## Matrix tracks

- `cooler` (`CoolerTrack`)
- `capcruncher` (`CapcruncherTrack`)
- `cooler_average` (`CoolerAverage`)

These require cooler-compatible input files.


## QuantNado tracks

- `quantnado_coverage` (`QuantNadoCoverageTrack`): single-sample coverage track.
- `quantnado_stranded_coverage` (`QuantNadoStrandedCoverageTrack`): mirrored +/− strand coverage.
- `quantnado_methylation` (`QuantNadoMethylationTrack`): per-CpG methylation scatter.
- `quantnado_variant` (`QuantNadoVariantTrack`): variant allele-frequency lollipop track.

These tracks accept either:

- a live `quantnado` object (or `dataset_path`) for region-time fetching, or
- precomputed xarray-like arrays with dims `(sample, position)`.
## Option reference

For full per-field options and descriptions:

- [Track Aliases](track_aliases.md)
- [Aesthetics Reference](aesthetics_reference.md)
- `GenomicFigure.track_options("<alias>")`
