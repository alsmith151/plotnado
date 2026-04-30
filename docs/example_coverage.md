# Example Coverage

This page maps track types to runnable examples and rendered outputs.

The default run of [examples/run_examples.py](https://github.com/alsmith151/plotnado/blob/main/examples/run_examples.py) replays the local example set. Remote Blueprint `bigwig_collection` / `bigwig_diff` coverage is opt-in with `python examples/run_examples.py --include-remote` because those scripts stream public BigWig files instead of checking large fixtures into the repo.

Rendered PNG outputs are written into `examples/output/` and synced into `docs/images/examples/`, so the gallery below shows the same results produced by the current code.

## Runnable script coverage

| Track / alias | Coverage level | Example path |
|---|---|---|
| `scalebar` | Full runnable | [examples/quickstart/01_first_plot.py](https://github.com/alsmith151/plotnado/blob/main/examples/quickstart/01_first_plot.py) |
| `axis` | Full runnable | [examples/quickstart/01_first_plot.py](https://github.com/alsmith151/plotnado/blob/main/examples/quickstart/01_first_plot.py) |
| `genes` | Full runnable | [examples/recipes/03_gene_label_strategies.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/03_gene_label_strategies.py) |
| `bigwig` | Full runnable | [examples/tracks/01_bigwig_styles.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/01_bigwig_styles.py) |
| `bed` | Full runnable | [examples/tracks/02_bed_and_narrowpeak.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/02_bed_and_narrowpeak.py) |
| `narrowpeak` | Full runnable | [examples/tracks/02_bed_and_narrowpeak.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/02_bed_and_narrowpeak.py) |
| `links` | Full runnable | [examples/tracks/03_links_annotations.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/03_links_annotations.py) |
| `hline` | Full runnable | [examples/tracks/03_links_annotations.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/03_links_annotations.py) |
| `vline` | Full runnable | [examples/tracks/03_links_annotations.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/03_links_annotations.py) |
| `bigwig_collection` | Runnable from remote Blueprint data | [examples/tracks/04_bigwig_collection_and_diff.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/04_bigwig_collection_and_diff.py) |
| `bigwig_diff` | Runnable from remote Blueprint data | [examples/tracks/04_bigwig_collection_and_diff.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/04_bigwig_collection_and_diff.py) |
| `cooler` | Full runnable | [examples/tracks/05_matrix_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/05_matrix_tracks.py) |
| `capcruncher` | Full runnable | [examples/tracks/05_matrix_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/05_matrix_tracks.py) |
| `cooler_average` | Full runnable | [examples/tracks/05_matrix_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/05_matrix_tracks.py) |
| `quantnado_coverage` | Full runnable | [examples/tracks/06_quantnado_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/06_quantnado_tracks.py) |
| `quantnado_stranded_coverage` | Full runnable | [examples/tracks/06_quantnado_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/06_quantnado_tracks.py) |
| `quantnado_methylation` | Full runnable | [examples/tracks/06_quantnado_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/06_quantnado_tracks.py) |
| `quantnado_variant` | Full runnable | [examples/tracks/06_quantnado_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/06_quantnado_tracks.py) |
| `highlight` | Full runnable | [examples/recipes/01_autoscale_overlay_highlight.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/01_autoscale_overlay_highlight.py) |
| `overlay` | Full runnable | [examples/recipes/01_autoscale_overlay_highlight.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/01_autoscale_overlay_highlight.py) |

## Runnable Example Gallery

### Basic figure

Minimal scale bar, axis, genes, and one in-memory signal track.

Source: [examples/basic_figure.py](https://github.com/alsmith151/plotnado/blob/main/examples/basic_figure.py)

![Basic figure](images/examples/basic_figure.png)

### Advanced features

Autocolor, highlight overlays, interval tracks, grouped autoscaling, and reference lines in one figure.

Source: [examples/advanced_features.py](https://github.com/alsmith151/plotnado/blob/main/examples/advanced_features.py)

![Advanced features](images/examples/advanced_features.png)

### Quickstart first plot

The minimal quickstart script used from the installation guide.

Source: [examples/quickstart/01_first_plot.py](https://github.com/alsmith151/plotnado/blob/main/examples/quickstart/01_first_plot.py)

![Quickstart first plot](images/examples/quickstart_first_plot.png)

### Alias and options discovery

Demonstrates `available_track_aliases()` and `track_options(...)` and saves a small multi-track plot.

Source: [examples/quickstart/02_aliases_and_options.py](https://github.com/alsmith151/plotnado/blob/main/examples/quickstart/02_aliases_and_options.py)

![Alias and options discovery](images/examples/quickstart_aliases_and_options.png)

### BigWig styles

Compares `fill`, `fragment`, `scatter`, and `std` BigWig render modes.

Source: [examples/tracks/01_bigwig_styles.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/01_bigwig_styles.py)

![BigWig styles](images/examples/track_bigwig_styles.png)

### BED and narrowPeak

Compares stacked BED intervals with score-colored narrowPeak rendering.

Source: [examples/tracks/02_bed_and_narrowpeak.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/02_bed_and_narrowpeak.py)

![BED and narrowPeak](images/examples/track_bed_and_narrowpeak.png)

### Links, hline, and vline

Shows BEDPE-style links together with horizontal and vertical annotation lines.

Source: [examples/tracks/03_links_annotations.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/03_links_annotations.py)

![Links, hline, and vline](images/examples/track_links_annotations.png)

### BigWig collection and diff

Fetches two public Blueprint plasma-cell RNA BigWigs once, stages tiny local crops for the target window in a temporary directory, then renders a `bigwig_collection` panel and a matching `bigwig_diff` panel from those local files.

Source: [examples/tracks/04_bigwig_collection_and_diff.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/04_bigwig_collection_and_diff.py)

![BigWig collection and diff](images/examples/track_bigwig_collection_diff.png)

### Cooler and cooler_average

Uses the checked-in GM12878 cooler fixture for a direct matrix view plus a simple averaged view built from the same file twice.

Source: [examples/tracks/05_matrix_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/05_matrix_tracks.py)

![Cooler and cooler_average](images/examples/track_matrix_cooler_average.png)

### CapCruncher viewpoint matrix

Uses the checked-in reporters-store fixture and resolves the `Slc25A37` viewpoint to the cooler group inside the HDF5 container.

Source: [examples/tracks/05_matrix_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/05_matrix_tracks.py)

![CapCruncher viewpoint matrix](images/examples/track_capcruncher.png)

### QuantNado coverage and stranded coverage

Uses local on-disk QuantNado fixtures to render real ATAC coverage and stranded RNA signal.

Source: [examples/tracks/06_quantnado_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/06_quantnado_tracks.py)

![QuantNado coverage and stranded coverage](images/examples/track_quantnado_signal.png)

### QuantNado methylation and variant tracks

The same script also renders file-backed methylation and variant panels from the compact test fixtures.

Source: [examples/tracks/06_quantnado_tracks.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/06_quantnado_tracks.py)

![QuantNado methylation](images/examples/track_quantnado_methylation.png)

![QuantNado variant](images/examples/track_quantnado_variant.png)

### Overlay, autoscale, and highlight

Demonstrates a shared overlay panel, grouped autoscaling, and region highlighting.

Source: [examples/recipes/01_autoscale_overlay_highlight.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/01_autoscale_overlay_highlight.py)

![Overlay, autoscale, and highlight](images/examples/recipe_autoscale_overlay_highlight.png)

### Theme, labels, and TOML round-trip

Applies the publication theme, saves a TOML figure definition, and reloads it to render again.

Source: [examples/recipes/02_theme_labels_toml.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/02_theme_labels_toml.py)

![Theme, labels, and TOML](images/examples/recipe_theme_labels.png)

Reloaded from TOML:

![Theme, labels, and TOML reloaded](images/examples/recipe_theme_labels_from_toml.png)

### Gene label strategies

Compares staggered, suppressed, and auto-expanded gene label placement strategies.

Source: [examples/recipes/03_gene_label_strategies.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/03_gene_label_strategies.py)

![Gene label strategies](images/examples/recipe_gene_label_strategies.png)

## Validate option coverage quickly

Use runtime metadata to inspect every field (track, aesthetics, label):

```python
from plotnado import GenomicFigure

for alias in sorted(GenomicFigure.available_track_aliases()):
    print(alias)
    GenomicFigure.track_options(alias)
```

For full generated tables, see [Aesthetics Reference](aesthetics_reference.md).
