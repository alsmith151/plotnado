# Track Catalog

This page summarizes the major track families, when to reach for them, and where to find runnable examples.

## Structural tracks

- `scalebar` (`ScaleBar`): genomic distance bar.
- `axis` (`GenomicAxis`): coordinate ticks and labels.
- `spacer` (`Spacer`): layout spacing.

## Signal tracks

- `bigwig` (`BigWigTrack`): continuous/interval signal (`style`: `fill`, `fragment`, `scatter`, `std`).
- `overlay` (`OverlayTrack`): multi-signal overlay panel with one shared y-axis derived from all component values.
- `bigwig_overlay` (`BigwigOverlay`): compatibility alias for `overlay`.
- `bigwig_collection` (`BigWigCollection`): multi-BigWig collection workflows.
- `bigwig_diff` (`BigWigDiff`): subtraction/ratio/log2ratio between two tracks.

Use `overlay` when the signals should be compared directly on one axis. Put `autoscale_group` on the overlay itself when it should match neighboring signal tracks, and use explicit `min_value` / `max_value` only when you want to pin that bound.

Scripts:

- [examples/tracks/01_bigwig_styles.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/01_bigwig_styles.py)
- [examples/recipes/01_autoscale_overlay_highlight.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/01_autoscale_overlay_highlight.py)

## Interval and annotation tracks

- `bed` (`BedTrack`): interval rectangles with optional labels.
- `narrowpeak` (`NarrowPeakTrack`): peak intervals with summit/score support.
- `genes` (`Genes`): gene models from bundled or custom annotations.
- `links` (`LinksTrack`): arcs for BEDPE-style interactions.
- `highlight` (`HighlightsFromFile`): region overlays.
- `hline` / `vline`: reference lines.

Scripts:

- [examples/tracks/02_bed_and_narrowpeak.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/02_bed_and_narrowpeak.py)
- [examples/tracks/03_links_annotations.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/03_links_annotations.py)

## Matrix tracks

- `cooler` (`CoolerTrack`)
- `capcruncher` (`CapcruncherTrack`)
- `cooler_average` (`CoolerAverage`)

These require cooler-compatible files.

## QuantNado tracks

- `quantnado_coverage`
- `quantnado_stranded_coverage`
- `quantnado_methylation`
- `quantnado_variant`

These support object-backed and array-backed workflows.

## Field-level options

- Runtime introspection: `GenomicFigure.track_options("<alias>")`
- Full generated table: [Aesthetics Reference](aesthetics_reference.md)
- Script mapping by track: [Example Coverage](example_coverage.md)

## Coverage snapshot

| Track / alias | Coverage | Primary example(s) |
|---|---|---|
| `scalebar`, `axis`, `genes` | Runnable script | [examples/quickstart/01_first_plot.py](https://github.com/alsmith151/plotnado/blob/main/examples/quickstart/01_first_plot.py), [examples/recipes/03_gene_label_strategies.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/03_gene_label_strategies.py) |
| `bigwig` | Runnable script | [examples/tracks/01_bigwig_styles.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/01_bigwig_styles.py), [examples/basic_figure.py](https://github.com/alsmith151/plotnado/blob/main/examples/basic_figure.py) |
| `bed`, `narrowpeak` | Runnable script | [examples/tracks/02_bed_and_narrowpeak.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/02_bed_and_narrowpeak.py) |
| `links`, `hline`, `vline` | Runnable script | [examples/tracks/03_links_annotations.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/03_links_annotations.py) |
| `highlight`, `overlay` | Runnable script | [examples/recipes/01_autoscale_overlay_highlight.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/01_autoscale_overlay_highlight.py) |
| `bigwig_collection`, `bigwig_diff` | Documented pattern | [Recipes](recipes.md) |
| `cooler`, `capcruncher`, `cooler_average` | Documented pattern | [Track Catalog](track_catalog.md) |
| QuantNado aliases | Documented pattern | [Track Construction](quickstart_tracks.md) |
