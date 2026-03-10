# Build Tracks Fast

Use aliases with `GenomicFigure.add_track()` or chain helper methods.

## Alias workflow

```python
from plotnado import GenomicFigure

fig = GenomicFigure()
fig.add_track("scalebar")
fig.add_track("axis")
fig.add_track("genes", genome="hg38")
fig.add_track("bigwig", data="signal.bw", title="ChIP", color="#1f77b4", alpha=0.7)
```

## Shorthand kwargs

You can pass many aesthetics and label options directly:

- Aesthetics fields (for example `color`, `alpha`, `style`) are routed to `aesthetics`.
- Label fields (for example `plot_title`, `title_location`) are routed to `label`.

## Discover aliases and options

```python
from plotnado import GenomicFigure

GenomicFigure.available_track_aliases()
GenomicFigure.track_options("bigwig")
GenomicFigure.track_options_markdown("genes")
```

Or from CLI:

```bash
plotnado track-options
plotnado track-options bigwig
plotnado track-options --all --output-format json
```

Runnable example: `python examples/quickstart/02_aliases_and_options.py`

Example output:

![Alias and options quickstart](images/examples/quickstart_aliases_and_options.png)
<<<<<<< HEAD

## QuantNado-backed tracks

Object-backed mode (fetch at plot time):

```python
from plotnado import GenomicFigure
from quantnado import QuantNado

qn = QuantNado.open("dataset.zarr")
fig = GenomicFigure()
fig.quantnado_coverage("sample1", quantnado=qn, title="Coverage")
fig.quantnado_methylation("sample1", quantnado=qn, title="Methylation")
fig.plot("chr1:100000-110000")
```

Array-backed mode (precomputed xarray-like arrays):

```python
fig = GenomicFigure()
fig.quantnado_variant(
    "sample1",
    allele_depth_ref_data=ref_da,
    allele_depth_alt_data=alt_da,
    title="Variant AF",
)
fig.plot("chr1:100000-110000")
```
=======
>>>>>>> 2273572867bad7bb6edf1bf8f5ecff6cd4752d5b
