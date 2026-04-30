# Aesthetics Reference

Use this page when you want the fastest route to track fields and styling options without reading source.

Every track exposes three groups of fields:

- `track`: constructor-level fields such as `title`, `height`, `autoscale_group`, and `color_group`
- `aesthetics`: visual controls such as `style`, `color`, `alpha`, `min_value`, and `max_value`
- `label`: title and scale annotation controls such as `plot_scale`, `title_location`, and `scale_location`

## Discover fields at runtime

```python
from plotnado import GenomicFigure, BigWigTrack, list_options

GenomicFigure.track_options("bigwig")
GenomicFigure.track_options("overlay")
GenomicFigure.track_options_markdown("genes")

BigWigTrack.options()
BigWigTrack.options_markdown()
list_options(BigWigTrack)
```

## High-value fields by category

| Category | Common fields | Notes |
| --- | --- | --- |
| Track | `title`, `height`, `autoscale_group`, `color_group` | These apply to the whole track panel |
| Aesthetics | `color`, `alpha`, `style`, `min_value`, `max_value` | These affect how the track is drawn |
| Label | `plot_title`, `plot_scale`, `title_location`, `scale_location` | These control title and scale text |

## Overlay-specific guidance

Overlay tracks are scaled as a single panel.

- `autoscale_group` belongs on the overlay itself when the whole panel should match sibling signal tracks.
- Explicit overlay `min_value` / `max_value` overrides grouped or global autoscale on that edge.
- `show_labels=False` is useful when the overlay is dense and the panel title carries enough meaning.

For behavior details, see [Aesthetics](aesthetics.md) and [Reference](reference.md).