# Track Aliases

`Figure.add_track()` accepts either:

- a concrete track instance (e.g. `BigWigTrack(data="a.bw")`), or
- a string alias (e.g. `"bigwig"`) plus constructor kwargs.

## Complete Alias Reference

| Alias | Track Class | Typical kwargs |
|---|---|---|
| `scalebar` | `ScaleBar` | `aesthetics={...}` |
| `scale` | `ScaleBar` | synonym of `scalebar` |
| `axis` | `GenomicAxis` | `aesthetics={...}` |
| `genes` | `Genes` | `genome="hg38"` or `data="genes.gtf.gz"` |
| `spacer` | `Spacer` | `height=...` |
| `bigwig` | `BigWigTrack` | `data="signal.bw"`, `aesthetics={"style": "fill"}` |
| `bed` | `BedTrack` | `data="regions.bed"` or in-memory `DataFrame` |
| `highlight` | `HighlightsFromFile` | `data="regions.bed"` |
| `bigwig_overlay` | `BigwigOverlay` | `tracks=["a.bw", "b.bw"]` |
| `bigwig_collection` | `BigWigCollection` | `files=[...]`, `aesthetics={"style": "overlay"}` |
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
