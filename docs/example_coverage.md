# Example Coverage

This page maps track types to runnable examples and rendered outputs.

All runnable scripts in [examples/run_examples.py](https://github.com/alsmith151/plotnado/blob/main/examples/run_examples.py) were re-run against the current codebase. Their PNG outputs are written into `examples/output/` and synced into `docs/images/examples/`, so the gallery below shows the same rendered results produced by the current code.

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

## Additional documented snippets

These tracks are documented with snippets and option discovery, but do not yet have dedicated runnable scripts in `examples/tracks/`.

### `bigwig_collection` / `bigwig_diff`

```python
gf.bigwig_collection(files=["a.bw", "b.bw"], title="Replicates")
gf.bigwig_diff(file_a="treated.bw", file_b="control.bw", title="Treated - Control")
```

### Matrix tracks

```python
gf.cooler(file="sample.mcool", resolution=10_000)
# gf.capcruncher(...)
# gf.cooler_average(...)
```

### QuantNado tracks

```python
gf.quantnado_coverage("sample1", quantnado=qn, title="Coverage")
gf.quantnado_methylation("sample1", quantnado=qn, title="Methylation")
```

## Validate option coverage quickly

Use runtime metadata to inspect every field (track, aesthetics, label):

```python
from plotnado import GenomicFigure

for alias in sorted(GenomicFigure.available_track_aliases()):
    print(alias)
    GenomicFigure.track_options(alias)
```

For full generated tables, see [Aesthetics Reference](aesthetics_reference.md).
