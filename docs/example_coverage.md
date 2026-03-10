# Example Coverage

This page maps track types to runnable examples and quick snippets.

## Runnable script coverage

| Track / alias | Coverage level | Example path |
|---|---|---|
| `scalebar` | Full runnable | `examples/quickstart/01_first_plot.py` |
| `axis` | Full runnable | `examples/quickstart/01_first_plot.py` |
| `genes` | Full runnable | `examples/recipes/03_gene_label_strategies.py` |
| `bigwig` | Full runnable | `examples/tracks/01_bigwig_styles.py` |
| `bed` | Full runnable | `examples/tracks/02_bed_and_narrowpeak.py` |
| `narrowpeak` | Full runnable | `examples/tracks/02_bed_and_narrowpeak.py` |
| `links` | Full runnable | `examples/tracks/03_links_annotations.py` |
| `hline` | Full runnable | `examples/tracks/03_links_annotations.py` |
| `vline` | Full runnable | `examples/tracks/03_links_annotations.py` |
| `highlight` | Full runnable | `examples/recipes/01_autoscale_overlay_highlight.py` |
| `overlay` | Full runnable | `examples/recipes/01_autoscale_overlay_highlight.py` |

## Additional documented snippets

These tracks are documented with snippets and option discovery, but do not yet have dedicated runnable scripts in `examples/tracks/`.

### `bigwig_collection` / `bigwig_diff`

```python
gf.bigwig_collection(files=["a.bw", "b.bw"], title="Replicates")
gf.bigwig_diff(file_a="treated.bw", file_b="control.bw", title="Treated - Control")
```

### Matrix tracks

```python
gf.cooler(file="sample.mcool", resolution=10_000)
# gf.capcruncher(...)
# gf.cooler_average(...)
```

### QuantNado tracks

```python
gf.quantnado_coverage("sample1", quantnado=qn, title="Coverage")
gf.quantnado_methylation("sample1", quantnado=qn, title="Methylation")
```

## Validate option coverage quickly

Use runtime metadata to inspect every field (track, aesthetics, label):

```python
from plotnado import GenomicFigure

for alias in sorted(GenomicFigure.available_track_aliases()):
    print(alias)
    GenomicFigure.track_options(alias)
```

For full generated tables, see [Aesthetics Reference](aesthetics_reference.md).
