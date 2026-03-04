# PlotNado

[![Tests](https://github.com/alsmith151/plotnado/actions/workflows/run_tests.yml/badge.svg)](https://github.com/alsmith151/plotnado/actions/workflows/run_tests.yml)

PlotNado is a lightweight Python package for building genome browser-style figures with a modern, chainable API.

> [!NOTE]
> This version of PlotNado is independent and does not require CoolBox.

## Install

```bash
pip install plotnado
```

For development from source:

```bash
git clone https://github.com/alsmith151/plotnado
cd plotnado
pip install -e .[dev,docs]
```

## Quick Start

```python
from plotnado import Figure
import numpy as np
import pandas as pd

bins = np.arange(1_000_000, 1_100_000, 1_000)
signal = pd.DataFrame({
    "chrom": "chr1",
    "start": bins,
    "end": bins + 1_000,
    "value": 5 + 2 * np.sin(np.linspace(0, 6, len(bins))),
})

fig = Figure()
fig.scalebar()
fig.axis()
fig.genes("hg38")
fig.bigwig(signal, title="Synthetic signal", style="fill", color="#1f77b4")
fig.save("quickstart.png", "chr1:1,010,000-1,080,000")
```

## Examples

- Start with `python examples/basic_figure.py`
- Then run `python examples/advanced_features.py`
- Full curated suite: `python examples/run_examples.py`
- Additional focused examples live in:
  - `examples/quickstart/`
  - `examples/tracks/`
  - `examples/recipes/`

All scripts write outputs to `examples/output/`.

## Key Features

- Chainable `Figure` API for fast composition.
- Alias-based track creation (`fig.add_track("bigwig", ...)`).
- Track option introspection at runtime.
- Built-in themes, autocolor, autoscale, and highlight overlays.
- Broad track support: BigWig, BED, narrowPeak, genes, axis, scalebar, links, `OverlayTrack`, cooler-based matrix tracks.

## Discover Track Options

```python
from plotnado import Figure, BigWigTrack

Figure.available_track_aliases()
Figure.track_options("bigwig")
Figure.track_options_markdown("bigwig")
BigWigTrack.options_markdown()
```

## Documentation

- Docs home: `docs/index.md`
- Quick start: `docs/quickstart.md`
- Track catalog + options: `docs/track_catalog.md`
- API and generated references: `docs/api_reference.md`

Build docs locally:

```bash
pip install -e .[docs]
mkdocs build --strict
mkdocs serve
```

## CLI

```python
plotnado plot chr1:1,000,000-1,100,000 -o browser_view.png
plotnado track-options bigwig
plotnado track-options --all --output-format json
```

