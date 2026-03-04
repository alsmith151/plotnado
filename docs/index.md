# PlotNado

PlotNado is a Python package for quick genome-browser-style plots with a Pydantic-first track architecture and a lightweight plotting API.

# Installation

To install PlotNado, you can use pip:

```bash
pip install git+https://github.com/alsmith151/plotnado

# Once this is on pypi, you will be able to install it with:
#pip install plotnado

```

# Usage

* Basic plotting [example notebook](basic_example.ipynb).
* Simple overlays of multiple tracks [example notebook](simple_overlays.ipynb).

## Track aliases

When using `Figure.add_track()`, you can pass aliases instead of constructing track classes directly.

| Alias | Track Class |
|---|---|
| `scalebar` / `scale` | `ScaleBar` |
| `axis` | `GenomicAxis` |
| `genes` | `Genes` |
| `spacer` | `Spacer` |
| `bigwig` | `BigWigTrack` |
| `bed` | `BedTrack` |
| `highlight` | `HighlightsFromFile` |
| `bigwig_overlay` | `BigwigOverlay` |
| `bigwig_collection` | `BigWigCollection` |
| `bigwig_diff` | `BigWigDiff` |
| `narrowpeak` | `NarrowPeakTrack` |
| `links` | `LinksTrack` |
| `hline` | `HLineTrack` |
| `vline` | `VLineTrack` |
| `cooler` | `CoolerTrack` |
| `capcruncher` | `CapcruncherTrack` |
| `cooler_average` | `CoolerAverage` |

For full details (constructor arguments, examples, and old→new migration mapping), see [Track Aliases](track_aliases.md).

## Shorthand Composition

`Figure` methods support shorthand composition for readability.

- Pass track fields directly (for example `title=...`, `height=...`).
- Pass aesthetics fields directly (for example `color=...`, `alpha=...`) and they are auto-packed into `aesthetics`.
- Pass label fields directly (for example `plot_title=False`) and they are auto-packed into `label`.

```python
from plotnado import Figure

fig = Figure()
fig.bigwig("chip.bw", title="ChIP", color="#1f77b4", alpha=0.7, plot_title=False)
```

Explicit nested models still work and can be mixed with shorthand; shorthand values take precedence for the keys provided.


