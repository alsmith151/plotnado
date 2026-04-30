# Recipes

Practical combinations of options and track types, with rendered output so you can see what the recipe actually produces.

## Autoscale + overlay + highlight

- Enable global autoscaling with `gf.autoscale(True)`.
- Put `autoscale_group` on the overlay itself when the overlay should sync with nearby signal panels.
- Mark regions with `gf.highlight(...)` and `gf.highlight_style(...)`.

![Autoscale + overlay + highlight](images/examples/recipe_autoscale_overlay_highlight.png)

Use this when you want one overlay panel to sit beside ordinary signal tracks without silently drifting onto a different y-scale.

Source: [examples/recipes/01_autoscale_overlay_highlight.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/01_autoscale_overlay_highlight.py)

## Theme + labels + TOML round-trip

- Apply publication defaults with `theme="publication"`.
- Control labels with shorthand kwargs (`title`, `title_color`, `title_location`).
- Save and reload figure definitions with TOML.

![Theme, labels, and TOML](images/examples/recipe_theme_labels.png)

Source: [examples/recipes/02_theme_labels_toml.py](https://github.com/alsmith151/plotnado/blob/main/examples/recipes/02_theme_labels_toml.py)

## Autocolor + color groups

- Enable auto palette assignment with `gf.autocolor()`.
- Use `color_group` so related tracks share one color.

```python
gf = GenomicFigure(theme="publication")
gf.autocolor()
gf.bed("THP1_peaks.bigBed", title="THP1 peaks", color_group="THP1")
gf.bigwig("THP1_signal.bw", title="THP1 signal", color_group="THP1")
gf.bed("K562_peaks.bigBed", title="K562 peaks", color_group="K562")
gf.bigwig("K562_signal.bw", title="K562 signal", color_group="K562")
```

## Track-style comparisons

- Compare BigWig styles (`fill`, `fragment`, `scatter`, `std`).
- Compare BED and narrowPeak interval behavior.

![BigWig style comparison](images/examples/track_bigwig_styles.png)

![BED and narrowPeak comparison](images/examples/track_bed_and_narrowpeak.png)

Scripts:

- [examples/tracks/01_bigwig_styles.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/01_bigwig_styles.py)
- [examples/tracks/02_bed_and_narrowpeak.py](https://github.com/alsmith151/plotnado/blob/main/examples/tracks/02_bed_and_narrowpeak.py)

More rendered runnable outputs, including quickstart, links/hline/vline, and gene-label examples, are collected on [Example Coverage](example_coverage.md).
