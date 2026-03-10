# Recipes

Practical combinations of options and track types.

## Autoscale + overlay + highlight

- Enable global autoscaling with `gf.autoscale(True)`.
- Overlay multiple signals with `gf.overlay(...)` / `OverlayTrack`.
- Mark regions with `gf.highlight(...)` and `gf.highlight_style(...)`.

Script: `examples/recipes/01_autoscale_overlay_highlight.py`

## Theme + labels + TOML round-trip

- Apply publication defaults with `theme="publication"`.
- Control labels with shorthand kwargs (`title`, `title_color`, `title_location`).
- Save and reload figure definitions with TOML.

Script: `examples/recipes/02_theme_labels_toml.py`

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

Scripts:

- `examples/tracks/01_bigwig_styles.py`
- `examples/tracks/02_bed_and_narrowpeak.py`
