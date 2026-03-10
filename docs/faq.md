# FAQ

## Do I need CoolBox?

No. PlotNado is independent and does not require CoolBox.

## Can I use DataFrames instead of files?

Yes. BigWig-like signals, BED-like intervals, highlights, and links can be created from in-memory DataFrames.

## How do I discover all options for a track?

Use runtime introspection:

```python
from plotnado import GenomicFigure
GenomicFigure.track_options("bigwig")
GenomicFigure.track_options_markdown("bigwig")
```

or CLI:

```bash
plotnado track-options bigwig
```

## How do I apply consistent style to all tracks?

`GenomicFigure()` uses publication defaults automatically. Use:

- `GenomicFigure(theme=...)` with `Theme.default()`, `Theme.minimal()`, or `Theme.publication()`.
- `GenomicFigure(theme=None)` to opt out of themed defaults.

## How do I share figure definitions?

Use TOML:

```python
fig.to_toml("plot.toml")
loaded = GenomicFigure.from_toml("plot.toml")
```
