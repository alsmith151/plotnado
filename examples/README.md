# PlotNado examples

This folder contains runnable, current-API examples.

## Start here

- `python examples/basic_figure.py`
- `python examples/advanced_features.py`
- `python examples/run_examples.py` runs the local example set.
- `python examples/run_examples.py --include-remote` also runs the Blueprint remote-BigWig example.

## Structured example sets

- `examples/quickstart/01_first_plot.py` - first figure in minutes.
- `examples/quickstart/02_aliases_and_options.py` - discover aliases/options and compose with shorthand.
- `examples/tracks/01_bigwig_styles.py` - fill/fragment/scatter BigWig rendering.
- `examples/tracks/02_bed_and_narrowpeak.py` - BED and narrowPeak interval tracks.
- `examples/tracks/03_links_annotations.py` - links plus horizontal/vertical annotation lines.
- `examples/tracks/04_bigwig_collection_and_diff.py` - remote Blueprint `bigwig_collection` and `bigwig_diff`.
- `examples/tracks/05_matrix_tracks.py` - local `cooler`, `cooler_average`, and `capcruncher` tracks.
- `examples/tracks/06_quantnado_tracks.py` - local QuantNado coverage, stranded, methylation, and variant tracks.
- `examples/recipes/01_autoscale_overlay_highlight.py` - autoscale, overlay, highlight styles.
- `examples/recipes/02_theme_labels_toml.py` - theme, label controls, TOML round-trip.
- `examples/recipes/03_gene_label_strategies.py` - compare gene label overlap strategies (stagger/suppress/auto-expand).

All scripts write output images into `examples/output/` by default.
