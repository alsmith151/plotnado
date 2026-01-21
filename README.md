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

For more examples, see the `examples/` directory.

