# FAQ

## Do I need CoolBox?

No. PlotNado is independent.

## Can I use DataFrames instead of files?

Yes. Signal and interval-style tracks support in-memory tabular inputs where applicable.

## How do I discover all options for a track?

```python
from plotnado import GenomicFigure

GenomicFigure.track_options("bigwig")
GenomicFigure.track_options_markdown("bigwig")
```

or CLI:

```bash
plotnado track-options bigwig
```

## What is the recommended coding style?

Use a single `GenomicFigure` instance (`gf`) and chain helper methods. See [Track Construction](quickstart_tracks.md).

## How do I share figure definitions?

```python
gf.to_toml("plot.toml")
loaded = GenomicFigure.from_toml("plot.toml")
```
