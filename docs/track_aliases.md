# Track Aliases

`Figure.add_track()` accepts either a concrete track instance or a string alias with kwargs.

```python
from plotnado import Figure

fig = Figure()
fig.add_track("scalebar", position="right")
fig.add_track("axis")
fig.add_track("genes", genome="hg38")
fig.add_track("bigwig", data="signal.bw", color="#1f77b4", alpha=0.7)
```

## Shorthand kwargs routing

- Track fields go directly to the track model.
- Aesthetics fields can be passed directly (auto-packed into `aesthetics`).
- Label fields can be passed directly (auto-packed into `label`).
- If both shorthand and nested models are provided, shorthand values win for overlapping keys.

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
| `overlay` | `OverlayTrack` | `tracks=[track_a, track_b]` or `tracks=["a.bw", "b.bw"]` |
| `bigwig_overlay` | `BigwigOverlay` | Backward-compatible alias of `overlay` |
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
fig.axis()
fig.genes("hg38")
fig.add_track("bigwig", data="chipseq.bw", title="ChIP", style="fill")
fig.add_track("vline", x_position=1_050_000)

fig.plot("chr1:1,000,000-1,100,000")
```

Use [Track Catalog](track_catalog.md) for track-specific examples and [API Reference](api_reference.md) for full option metadata.
