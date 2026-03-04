# Recipes

Practical combinations of options and track types.

## Autoscale + overlay + highlight

- Enable global autoscaling with `fig.autoscale(True)`.
- Overlay multiple signals with `fig.overlay(...)` / `OverlayTrack`.
- Mark regions with `fig.highlight(...)` and `fig.highlight_style(...)`.

Script: `examples/recipes/01_autoscale_overlay_highlight.py`

![Autoscale overlay highlight](images/examples/recipe_autoscale_overlay_highlight.png)

## Theme + labels + TOML round-trip

- Apply publication defaults with `Theme.publication()`.
- Control label rendering with shorthand label kwargs.
- Save and reload figure definitions with TOML.

Script: `examples/recipes/02_theme_labels_toml.py`

![Theme labels recipe](images/examples/recipe_theme_labels.png)
![Theme labels from TOML](images/examples/recipe_theme_labels_from_toml.png)

## Track-style comparisons

- Compare BigWig styles (`fill`, `fragment`, `scatter`).
- Compare BED and narrowPeak interval behavior.

Scripts:

- `examples/tracks/01_bigwig_styles.py`
- `examples/tracks/02_bed_and_narrowpeak.py`

![Track style comparison](images/examples/track_bigwig_styles.png)
