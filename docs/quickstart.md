# Quick Start

This guide gets you from installation to a first saved figure quickly.

PlotNado has two common entry points:

- Use the Python API when you want figure construction in code.
- Use the CLI when you want an editable YAML template and a file-driven workflow.

## 1) Build a minimal figure

```python
from plotnado import GenomicFigure
import numpy as np
import pandas as pd

bins = np.arange(1_000_000, 1_100_000, 1_000)
signal = pd.DataFrame(
    {
        "chrom": "chr1",
        "start": bins,
        "end": bins + 1_000,
        "value": 5 + 2 * np.sin(np.linspace(0, 6, len(bins))),
    }
)

gf = GenomicFigure(theme="publication")
gf.scalebar()
gf.axis()
gf.genes("hg38")
gf.bigwig(signal, title="Synthetic signal", style="fill")
gf.save("quickstart.png", "chr1:1,010,000-1,080,000")
```

## 2) Run a ready-made script

```bash
python examples/quickstart/01_first_plot.py
```

## 3) Generate a template with the CLI

```bash
plotnado init sample1.bw sample2.bw peaks.narrowpeak --auto --output template.yaml
plotnado validate template.yaml
plotnado plot template.yaml --region chr1:1,000,000-1,100,000 --output quickstart-cli.png
```

`plotnado plot` also accepts gene symbols such as `--region GNAQ` when the template has a `genome` value.

If you want to stay in Python after generating a template:

```python
from plotnado import GenomicFigure

gf = GenomicFigure.from_template("template.yaml")
gf.save("quickstart-from-template.png", region="chr1:1,000,000-1,100,000")
```

## 4) Start from a public UCSC hub in a notebook

PlotNado can flatten supported tracks from a UCSC-style hub directly into a figure. In notebooks, `track_visibility_widget()` is the fastest way to hide or re-enable tracks without rebuilding the figure.

```python
from plotnado import GenomicFigure

hub = "https://ftp.ebi.ac.uk/pub/databases/blueprint/releases/current_release/homo_sapiens/hub/hub.txt"
gf = GenomicFigure.from_ucsc_hub(hub, genome="hg38", include_hidden=False)

widget = gf.track_visibility_widget("chr12:6908053-6997143")
widget
```

Install the notebook extra first if you want the widget controls:

```bash
uv pip install 'plotnado[notebook]'
```

Set `include_hidden=True` when you want hub tracks that are hidden by default to be preserved as zero-height tracks that can be turned back on later.

## 5) Next steps

- Learn all track-add patterns: [Track Construction](quickstart_tracks.md)
- Learn the template workflow: [CLI](cli.md)
- Review practical defaults: [Best Practices](best_practices.md)
- Explore track families: [Track Catalog](track_catalog.md)
- See a longer live notebook walkthrough: [Worked Examples](worked_examples.ipynb)
