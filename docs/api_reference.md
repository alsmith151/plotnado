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

## `GenomicFigure` helper methods: automatic kwargs

For helper methods like `gf.bigwig(...)`, kwargs can be provided in shorthand form and PlotNado routes them automatically:

- Track fields: passed directly.
- Aesthetics fields: routed into `aesthetics`.
- Label fields: routed into `label`.

```python
from plotnado import GenomicFigure

gf = GenomicFigure()
gf.bigwig(
    "signal.bw",
    title="Sample A",      # label
    title_color="black",   # label
    style="std",           # aesthetics
    color="#1f77b4",       # aesthetics
    alpha=0.8,             # aesthetics
)
```

`color_group` is a track-level kwarg and works well with `gf.autocolor()` for consistent sample coloring:

```python
gf = GenomicFigure(theme="publication")
gf.autocolor()
gf.bed("sampleA.bigBed", title="A peaks", color_group="sampleA")
gf.bigwig("sampleA.bw", title="A signal", color_group="sampleA")
```

## Common entry points

- `plotnado.GenomicFigure`: high-level composition (`add_track`, `plot`, `plot_regions`, `plot_gene`, `from_template`, `from_igv_session`, `to_toml`, `from_toml`).
- `plotnado.Template`: YAML model used by the CLI and `GenomicFigure.from_template()`.
- `plotnado.TemplateCompiler`: converts a `Template` into a reusable render plan.
- `plotnado.Theme`: built-in or custom visual defaults.
- `plotnado.tracks.*`: concrete track classes when you want explicit model construction.
- `plotnado.parse_igv_session`: parse an IGV session XML into an `IgvSession` (`.template`, `.locus`, `.genome`).

## Track editing

After a figure is built, tracks can be edited in-place.

| Method | Description |
| --- | --- |
| `fig["title"]` / `fig[i]` | Access a track by title (case-insensitive) or index |
| `fig.update_track(key, **kw)` | Update fields on one track |
| `fig.update_track(**kw)` | Update all tracks (no key) |
| `fig.update_track(track_type=…, **kw)` | Bulk update filtered by type |
| `fig.update_track(group=…, **kw)` | Bulk update filtered by autoscale group |
| `fig.update_track(where=fn, **kw)` | Bulk update with a predicate |
| `fig.remove_track(key)` | Remove a track by title or index |
| `fig.add_track(…, position="top")` | Prepend instead of append |

All editing methods return `self` for chaining.

For practical usage, prefer [Quick Start](quickstart.md), [Track Catalog](track_catalog.md), and [Recipes](recipes.md).
