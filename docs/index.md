# PlotNado

PlotNado is a Python package for making clean, publication-ready genome browser plots with a fast, chainable API.

## Install

```bash
pip install plotnado
```

## 5-minute quick start

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

fig = GenomicFigure()
fig.scalebar().axis().genes("hg38")
fig.bigwig(signal, title="Synthetic signal", style="fill")
fig.save("quickstart.png", "chr1:1,010,000-1,080,000")
```

Example output:

![Quickstart plot](images/examples/quickstart_first_plot.png)

## Learn by task

- New user setup: [Installation](installation.md)
- First figure quickly: [Quick Start](quickstart.md)
- Add tracks faster with aliases: [Build Tracks Fast](quickstart_tracks.md)
- Browse track types and options: [Track Catalog](track_catalog.md)
- Multi-panel plotting and figure workflows: [Figure Workflows](figure_workflows.md)
- Practical configurations: [Recipes](recipes.md)
- Input data formats and requirements: [Data Inputs](data_inputs.md)
- CLI usage and option discovery: [CLI](cli.md)

## Example scripts

Run these from the repository root:

```bash
python examples/basic_figure.py
python examples/advanced_features.py
python examples/run_examples.py
```

Focused examples are organized in:

- `examples/quickstart/`
- `examples/tracks/`
- `examples/recipes/`

Sample outputs:

![Basic figure](images/examples/basic_figure.png)
![Advanced features](images/examples/advanced_features.png)

## Option discovery

```python
from plotnado import GenomicFigure

GenomicFigure.available_track_aliases()
GenomicFigure.track_options("bigwig")
GenomicFigure.track_options_markdown("genes")
```


