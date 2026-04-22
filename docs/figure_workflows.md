# Figure Workflows

## Single region plot

```python
gf.plot("chr1:1,000,000-1,100,000")
```

## Extend around region

```python
gf.plot("chr1:1,000,000-1,050,000", extend=0.25)
```

## Plot many regions

```python
gf.plot_regions(["chr1:1,000,000-1,050,000", "chr1:1,060,000-1,110,000"], ncols=2)
```

`plot_regions` also accepts a BED path.

## Plot by gene symbol

```python
gf.plot_gene("GNAQ", extend=0.5)
```

## Load from an IGV session

Convert a saved IGV `.xml` session directly into a `GenomicFigure`:

```python
from plotnado import GenomicFigure

fig, locus = GenomicFigure.from_igv_session("session.xml")
fig.plot(locus)  # locus is the session's saved view region
```

Track colors, fixed scales (`min`/`max`), autoscale groups, and gene annotations are
all preserved. The session locus is returned alongside the figure so you can use it
directly or override it.

For lower-level access, `parse_igv_session` returns the intermediate `IgvSession`
object (`.template`, `.locus`, `.genome`):

```python
from plotnado import parse_igv_session, GenomicFigure

session = parse_igv_session("session.xml")
fig = GenomicFigure.from_template(session.template)
fig.plot(session.locus)
```

## Edit tracks after building a figure

Access and mutate individual tracks by title (case-insensitive) or index:

```python
fig["25_NK"].color = "steelblue"
fig[0].height = 2.0
```

Update a single track by key:

```python
fig.update_track("CAT-MV411_H3K27ac", color="purple", max_value=500)
```

Update multiple tracks at once — omit `key` and optionally filter:

```python
fig.update_track(height=0.2)                                    # all tracks
fig.update_track(track_type="bigwig", height=0.3)               # by type
fig.update_track(group="igv_group_3", max_value=1000.0)         # by autoscale group
fig.update_track(where=lambda t: t.title and "CD" in t.title, color="coral")
```

Remove a track:

```python
fig.remove_track("25_NK")
fig.remove_track(0)
```

Add a track at the top or bottom (default):

```python
fig.bigwig("new.bw", title="extra", position="top")
fig.genes("hg38", position="top")
```

All editing methods return `self` for chaining:

```python
(fig
    .remove_track("01_HSC")
    .update_track(track_type="bigwig", height=0.25)
    .bigwig("new.bw", title="extra", position="top")
)
```

## Save and load figure configs

```python
gf.to_toml("figure.toml")
loaded = GenomicFigure.from_toml("figure.toml")
loaded.save("figure.png", "chr1:1,000,000-1,100,000")
```

See `examples/recipes/02_theme_labels_toml.py` for a full round-trip script.
