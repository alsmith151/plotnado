# API Reference

This page links to generated API reference content.

- [Reference Summary](reference/SUMMARY.md)

## Runtime Option Discovery

For notebook and REPL usage, track options are introspectable at runtime:

```python
from plotnado import BigWigTrack, list_options

# Class-level helper
BigWigTrack.options()
BigWigTrack.options_markdown()

# Module-level helper
list_options(BigWigTrack)
```

Both return a dictionary split into:
- `track`: core track fields
- `aesthetics`: flattened aesthetic fields accepted directly in constructors
- `label`: unified label controls (`track.label` / `label={...}`)

For alias-based workflows:

```python
from plotnado import Figure
from plotnado import Theme

Figure.available_track_aliases()
Figure.track_options("bigwig")
Figure.track_options_markdown("bigwig")

theme = Theme.publication()
fig = Figure(theme=theme)
```

Current aesthetics coverage (derived from `Track.options()`):

For the always up-to-date generated table per track (track/aesthetics/label fields), see [Aesthetics Reference](aesthetics_reference.md).

- `BigWigTrack`: `alpha`, `color`, `fill`, `linewidth`, `max_value`, `min_value`, `scatter_point_size`, `style`
- `BedTrack`: `alpha`, `color`, `display`, `edge_color`, `font_size`, `interval_height`, `label_field`, `max_rows`, `rect_linewidth`, `show_labels`
- `Genes`: `alpha`, `arrow_size`, `color`, `display`, `exon_edge_color`, `exon_linewidth`, `fill`, `gene_label_font_size`, `gene_label_style`, `interval_height`, `intron_color`, `intron_linewidth`, `max_number_of_rows`, `minimum_gene_length`, `style`
- `ScaleBar`: `bar_linewidth`, `color`, `font_size`, `label_offset`, `position`, `scale_distance`, `style`, `tick_height`, `tick_linewidth`
- `GenomicAxis`: `axis_linewidth`, `chromosome_fontweight`, `color`, `font_size`, `num_ticks`, `show_chromosome`, `tick_color`, `tick_height`, `tick_linewidth`
- `HighlightsFromFile`: `alpha`, `color`, `edge_color`, `linewidth`
- `LinksTrack`: `alpha`, `cmap`, `color`, `color_by_score`, `edge_color`, `linewidth`, `max_height`, `max_score`, `min_score`, `y_baseline`
- `BigwigOverlay`: `alpha`, `colors`, `max_value`, `min_value`, `show_labels`
- `BigWigCollection`: `alpha`, `colors`, `labels`, `style`
- `BigWigDiff`: `bar_alpha`, `linewidth`, `negative_color`, `positive_color`, `zero_line_alpha`, `zero_line_color`, `zero_line_width`
- `CoolerTrack`: `cmap`, `max_value`, `min_value` (legacy aliases `vmax`/`vmin` still supported)
- `NarrowPeakTrack`: `alpha`, `cmap`, `color`, `color_by`, `display`, `edge_color`, `font_size`, `interval_height`, `label_field`, `max_rows`, `max_score`, `min_score`, `rect_linewidth`, `show_labels`, `show_summit`, `summit_color`, `summit_width`

If the summary is empty, ensure the reference generator script is configured for the desired module paths.
