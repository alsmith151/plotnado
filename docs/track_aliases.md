# Track Aliases

`GenomicFigure.add_track()` accepts either a concrete track object or a string alias.

## Alias usage

```python
from plotnado import GenomicFigure

gf = GenomicFigure()
gf.add_track("scalebar")
gf.add_track("axis")
gf.add_track("genes", genome="hg38")
gf.add_track("bigwig", data="signal.bw", style="fill")
```

## How kwargs are routed

- Track constructor fields: passed directly.
- Aesthetics fields (for example `color`, `alpha`, `style`): routed to `aesthetics`.
- Label fields (for example `title`, `title_color`, `title_location`): routed to `label`.
- If both shorthand and nested objects are supplied, shorthand values win for overlapping keys.

## Common aliases

- Structural: `scalebar`, `axis`, `genes`, `spacer`
- Signal: `bigwig`, `overlay`, `bigwig_collection`, `bigwig_diff`
- Interval/annotation: `bed`, `narrowpeak`, `links`, `highlight`, `hline`, `vline`
- Matrix: `cooler`, `capcruncher`, `cooler_average`
- QuantNado: `quantnado_coverage`, `quantnado_stranded_coverage`, `quantnado_methylation`, `quantnado_variant`

## Get the full, current alias map

```python
from plotnado import GenomicFigure

aliases = GenomicFigure.available_track_aliases()
print("\n".join(f"{k} -> {v}" for k, v in sorted(aliases.items())))
```

Use [Track Catalog](track_catalog.md) for track families and [Aesthetics Reference](aesthetics_reference.md) for full field-level options.
