# Track Aliases

`Figure.add_track()` accepts either:

- a concrete track instance (e.g. `BigWigTrack(data="a.bw")`), or
- a string alias (e.g. `"bigwig"`) plus constructor kwargs.

## Shorthand Composition

Alias-based composition supports ergonomic kwargs routing:

- Track fields are passed through directly.
- Aesthetics fields are accepted directly and routed into `aesthetics`.
- Label fields are accepted directly and routed into `label`.

```python
from plotnado import Figure

fig = Figure()
fig.add_track("bigwig", data="a.bw", color="#1f77b4", alpha=0.7, plot_title=False)
fig.add_track("scalebar", position="right", scale_size=10)
```

You can still pass explicit nested models (for example `aesthetics=...`); if both are provided, shorthand keys override the matching nested keys.

## Complete Alias Reference

| Alias | Track Class | Typical kwargs |
|---|---|---|
| `scalebar` | `ScaleBar` | `position="left"`, `scale_distance=...` |
| `scale` | `ScaleBar` | synonym of `scalebar` |
| `axis` | `GenomicAxis` | `num_ticks=...`, `show_chromosome=True` |
| `genes` | `Genes` | `genome="hg38"` or `data="genes.gtf.gz"` |
| `spacer` | `Spacer` | `height=...` |
| `bigwig` | `BigWigTrack` | `data="signal.bw"`, `style="fill"`, `color=...` |
| `bed` | `BedTrack` | `data="regions.bed"` or in-memory `DataFrame` |
| `highlight` | `HighlightsFromFile` | `data="regions.bed"` |
| `bigwig_overlay` | `BigwigOverlay` | `tracks=["a.bw", "b.bw"]` |
| `bigwig_collection` | `BigWigCollection` | `files=[...]`, `style="overlay"` |
| `bigwig_diff` | `BigWigDiff` | `file_a="a.bw"`, `file_b="b.bw"`, `method="subtract"` |
| `narrowpeak` | `NarrowPeakTrack` | `data="peaks.narrowPeak"` |
| `links` | `LinksTrack` | `data="loops.bedpe"` |
| `hline` | `HLineTrack` | `y_value=...` |
| `vline` | `VLineTrack` | `x_position=...` |
| `cooler` | `CoolerTrack` | `file="matrix.mcool"`, `resolution=10000` |
| `capcruncher` | `CapcruncherTrack` | `file=...`, `viewpoint=...`, `normalisation=...` |
| `cooler_average` | `CoolerAverage` | `files=[... ]`, `resolution=...` |

## Example

```python
from plotnado import Figure

fig = Figure(width=12)
fig.add_track("axis")
fig.add_track("genes", genome="hg38")
fig.add_track("bigwig", data="chipseq.bw", title="ChIP")
fig.add_track("bigwig_diff", file_a="condA.bw", file_b="condB.bw", method="log2ratio")
fig.add_track("vline", x_position=1050000)

fig.plot("chr1:1,000,000-1,100,000")
```

## Migration Map (Old → New)

| Old (CoolBox-era) | New |
|---|---|
| `MatrixCapcruncher` | `CoolerTrack` |
| `MatrixCapcruncherAverage` | `CoolerAverage` |
| `BigwigFragment` | `BigWigTrack(style="fragment")` |
| `BigwigFragmentCollection` | `BigWigCollection` |
| `BigwigSubtraction` | `BigWigDiff` |
| `BedMemory` | `BedTrack(data=df)` |
| `Genes` / `GeneTrack` | `Genes` / `GeneTrack` |
| `GenomicAxis` | `AxisTrack` |
| `ScaleBar` | `ScaleBarTrack` |
| `Spacer` | `SpacerTrack` |
| `HighlightAnnotation` | `HighlightsFromFile` |
