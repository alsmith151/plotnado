# Figure Workflows

## Single region plot

```python
fig.plot("chr1:1,000,000-1,100,000")
```

## Extend around region

```python
fig.plot("chr1:1,000,000-1,050,000", extend=0.25)
```

## Plot many regions

```python
fig.plot_regions(["chr1:1,000,000-1,050,000", "chr1:1,060,000-1,110,000"], ncols=2)
```

`plot_regions` also accepts a BED path.

## Plot by gene symbol

```python
fig.plot_gene("DDX11L1", extend=0.5)
```

## Save and load figure configs

```python
fig.to_toml("figure.toml")
loaded = Figure.from_toml("figure.toml")
loaded.save("figure.png", "chr1:1,000,000-1,100,000")
```

See `examples/recipes/02_theme_labels_toml.py` for a full round-trip script.

Workflow-style output example:

![Advanced workflow figure](images/examples/advanced_features.png)
