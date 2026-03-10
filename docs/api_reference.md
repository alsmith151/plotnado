# API Reference

Use this section when you need API-level details beyond the task-oriented guides.

- Generated module/class reference: [Reference Summary](reference/SUMMARY.md)
- Auto-generated options reference for every track: [Aesthetics Reference](aesthetics_reference.md)

## Discover options at runtime

```python
from plotnado import GenomicFigure, BigWigTrack, list_options

GenomicFigure.available_track_aliases()
GenomicFigure.track_options("bigwig")
GenomicFigure.track_options_markdown("bigwig")

BigWigTrack.options()
BigWigTrack.options_markdown()
list_options(BigWigTrack)
```

Each option payload is split into:

- `track`: top-level constructor fields.
- `aesthetics`: style model fields (`aesthetics={...}` or shorthand kwargs).
- `label`: label controls (`label={...}` or shorthand kwargs).

## Common entry points

- `plotnado.GenomicFigure`: high-level composition (`add_track`, `plot`, `plot_regions`, `plot_gene`, `to_toml`, `from_toml`).
- `plotnado.Theme`: built-in or custom visual defaults.
- `plotnado.tracks.*`: concrete track classes when you want explicit model construction.

For practical usage, prefer [Quick Start](quickstart.md), [Track Catalog](track_catalog.md), and [Recipes](recipes.md).
