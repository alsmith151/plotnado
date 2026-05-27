# PlotNado

PlotNado creates genome browser-style figures from Python or YAML templates.
It is designed for reproducible genomic plots in scripts, notebooks, Quarto docs, and CLI workflows.
Use it when you want compact tracks for BigWig-like signals, BED/narrowPeak intervals, links, genes, and optional matrix/QuantNado data.

<p align="center">
  <img src="docs/images/Logo.jpeg" alt="PlotNado logo" width="220">
</p>

<p align="center">
  <a href="https://github.com/alsmith151/plotnado/actions/workflows/run_tests.yml">
    <img alt="Tests" src="https://github.com/alsmith151/plotnado/actions/workflows/run_tests.yml/badge.svg">
  </a>
</p>

## Start

Install for day-to-day use:

```bash
uv tool install plotnado
plotnado --help
```

Work from source:

```bash
git clone https://github.com/alsmith151/plotnado
cd plotnado
uv sync --extra dev --extra docs
uv run pytest tests/
```

## Python API

```python
from plotnado import GenomicFigure
import numpy as np
import pandas as pd

bins = np.arange(1_000_000, 1_100_000, 1_000)
signal = pd.DataFrame({
    "chrom": "chr1",
    "start": bins,
    "end": bins + 1_000,
    "value": 5 + 2 * np.sin(np.linspace(0, 6, len(bins))),
})

fig = (
    GenomicFigure()
    .scalebar()
    .axis()
    .bigwig(signal, title="Synthetic signal", style="fill", color="#1f77b4")
)
fig.save("quickstart.png", region="chr1:1,010,000-1,080,000")
```

## CLI + YAML

```bash
uv run plotnado init sample1.bw sample2.bw peaks.narrowpeak --auto --output template.yaml
uv run plotnado validate template.yaml
uv run plotnado plot template.yaml --region chr1:1,000,000-1,100,000 --output browser_view.png
```

Templates are editable YAML:

```yaml
genome: hg38
guides:
  genes: true
tracks:
  - path: sample1.bw
    type: bigwig
    title: sample1
  - path: peaks.narrowpeak
    type: narrowpeak
    title: peaks
```

The same template can be loaded from Python:

```python
from plotnado import GenomicFigure

fig = GenomicFigure.from_template("template.yaml")
fig.save("browser_view.png", region="chr1:1,000,000-1,100,000")
```

## What It Covers

- Chainable `GenomicFigure` API for programmatic plotting.
- `plotnado init`, `plotnado validate`, and `plotnado plot` for YAML workflows.
- In-memory tabular examples for reproducible docs and notebooks.
- Runtime alias and option lookup through `GenomicFigure.available_track_aliases()` and `GenomicFigure.track_options(...)`.
- Optional cooler, CapCruncher, and QuantNado tracks when their dependencies and input datasets are installed.

## Documentation

- [Installation](docs/installation.qmd)
- [Quick Start](docs/quickstart.qmd)
- [Track Catalog](docs/track_catalog.qmd)
- [Aesthetics](docs/aesthetics.qmd)
- [CLI](docs/cli.qmd)
- [Troubleshooting](docs/troubleshooting.qmd)

## Troubleshooting

If a plot is empty, first check that the region overlaps your data and that chromosome names match (`chr1` versus `1`).
If Quarto uses the wrong Python, render with `QUARTO_PYTHON=.venv/bin/python quarto render`.
For option names, use `GenomicFigure.track_options("bigwig")`.

## Development

```bash
uv sync --extra dev --extra docs
uv run pytest tests/
uv run python examples/run_examples.py
uv run plotnado --help
QUARTO_PYTHON=.venv/bin/python quarto render
```

Docs are Quarto pages with executable Python cells. Prefer deterministic in-memory data in docs examples and keep generated option tables secondary to plotted examples.

Contributions should start with a source install, the test suite, the example runner, and a Quarto render. The main branch is protected, so open a pull request for review.

## License

PlotNado is licensed under GPL-3.0-or-later. See [LICENSE](LICENSE).
