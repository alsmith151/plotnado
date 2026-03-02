# PlotNado

[![Tests](https://github.com/alsmith151/plotnado/actions/workflows/run_tests.yml/badge.svg)](https://github.com/alsmith151/plotnado/actions/workflows/run_tests.yml)

PlotNado is a lightweight Python package for creating genome browser-style plots with a focus on simplicity, rich aesthetics, and performance.

> [!NOTE]
> This new version of PlotNado is independent and does not require CoolBox.

## Key Features

- **Rich Aesthetics**: Premium default designs with support for custom palettes and styles.
- **Multiple Track Types**:
  - `BigWigTrack`: Signal data visualization with stairs or scatter styles.
  - `Genes`: Intelligent gene annotation display with intron/exon structure and collision-aware labels.
  - `BedTrack`: Simple and clean BED record display.
  - `GenomicAxis`: Clear x-axis with genomic coordinates.
  - `Highlights`: Highlight specific regions across all tracks.
- **Advanced Plotting**:
  - `BigwigOverlay`: Overlay multiple signals on a single axis.
  - `Autoscaler`: Automatically share y-axis scales across tracks.
  - `Figure API`: High-level API for easy composition, automatic coloring, and region manipulation.

## Installation

```bash
pip install git+https://github.com/alsmith151/plotnado
```

## Build Documentation

```bash
python -m pip install -r requirements-docs.txt
python -m mkdocs build --strict
```

## Basic Usage

```python
from plotnado import Figure

fig = Figure()
fig.add_track('scalebar')
fig.add_track('axis')
fig.add_track('genes', genome='hg38')
fig.add_track('bigwig', data='signal.bw', title='ChIP-seq')

fig.plot('chr1:1,000,000-1,100,000')
```

## Advanced Usage

```python
fig = Figure()
fig.autoscale(True).autocolor(palette='Set1')
fig.highlight('chr1:1045000-1055000')

fig.add_track('bigwig_overlay', tracks=[...])
fig.plot('chr1:1000000-1100000')
```

## Discoverability in IDE and notebooks

```python
from plotnado import BigWigTrack, Figure, LabelConfig

# Programmatic option discovery
BigWigTrack.options()            # dict with track/aesthetics/label sections
BigWigTrack.options_markdown()   # notebook-friendly markdown table

# Alias-based discovery
Figure.available_track_aliases()
Figure.track_options("bigwig")

# Unified label config
track = BigWigTrack(
  color="#1f77b4",
  label=LabelConfig(plot_title=True, plot_scale=True, label_box_alpha=0.8),
)
```

## Track Aliases

`Figure.add_track()` accepts either a track instance or a string alias. The following aliases are supported:

Full docs page: `docs/track_aliases.md`.

| Alias | Track Class | Notes |
|---|---|---|
| `scalebar` | `ScaleBar` | Scale bar track |
| `scale` | `ScaleBar` | Synonym for `scalebar` |
| `axis` | `GenomicAxis` | Genomic coordinate axis |
| `genes` | `Genes` | Gene annotations (`genome` or `data`) |
| `spacer` | `Spacer` | Empty spacing track |
| `bigwig` | `BigWigTrack` | BigWig signal track |
| `bed` | `BedTrack` | BED/BigBed intervals |
| `highlight` | `HighlightsFromFile` | Region highlighting from BED/BigBed |
| `bigwig_overlay` | `BigwigOverlay` | Overlay multiple BigWig tracks |
| `bigwig_collection` | `BigWigCollection` | Collection track (`overlay`/`stacked`) |
| `bigwig_diff` | `BigWigDiff` | Difference/ratio between two BigWig tracks |
| `narrowpeak` | `NarrowPeakTrack` | NarrowPeak intervals |
| `links` | `LinksTrack` | BEDPE-style interaction arcs |
| `hline` | `HLineTrack` | Horizontal annotation line |
| `vline` | `VLineTrack` | Vertical annotation line |
| `cooler` | `CoolerTrack` | Cooler matrix heatmap (`.cool`/`.mcool`) |
| `capcruncher` | `CapcruncherTrack` | CapCruncher-style Cooler track |
| `cooler_average` | `CoolerAverage` | Average matrix across multiple Cooler files |

Example:

```python
fig = Figure()
fig.add_track("cooler", file="sample.mcool", resolution=10000)
fig.add_track("bigwig_diff", file_a="a.bw", file_b="b.bw", method="log2ratio")
fig.add_track("vline", x_position=1050000)
```

For more examples, see the `examples/` directory.

