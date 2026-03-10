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

## Save and load figure configs

```python
gf.to_toml("figure.toml")
loaded = GenomicFigure.from_toml("figure.toml")
loaded.save("figure.png", "chr1:1,000,000-1,100,000")
```

See `examples/recipes/02_theme_labels_toml.py` for a full round-trip script.
