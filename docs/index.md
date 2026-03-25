# PlotNado

![PlotNado logo](assets/plotnado-logo.svg)

PlotNado is a Python library for clean, publication-ready genome browser figures.

It exposes two parallel interfaces:

- A fluent Python API for programmatic figure construction
- A CLI that turns genomic files into editable YAML templates and rendered plots

## Install

```bash
uv venv
source .venv/bin/activate
uv pip install plotnado
```

## Preferred workflow style

```python
from plotnado import GenomicFigure

gf = GenomicFigure(theme="publication")
gf.autocolor()
gf.scalebar()
gf.genes("hg38", display="expanded", minimum_gene_length=1e5)

gf.bigwig("signal_1.bw", title="H3K4me1", color_group="H3K4me1", style="std")
gf.bigwig("signal_2.bw", title="H3K4me3", color_group="H3K4me3", style="std")

gf.axis()
gf.plot_gene("GNAQ")
```

## Template-driven workflow

```bash
plotnado init *.bw peaks/*.narrowpeak --auto --output template.yaml
plotnado validate template.yaml
plotnado plot template.yaml --region chr1:1,000,000-1,100,000 --output region.png
```

Templates can also be consumed from Python with `GenomicFigure.from_template("template.yaml")`.

## Read by task

- New environment setup: [Installation](installation.md)
- First successful figure: [Quick Start](quickstart.md)
- Ways to add tracks: [Track Construction](quickstart_tracks.md)
- Production guidance: [Best Practices](best_practices.md)
- Track families and support matrix: [Track Catalog](track_catalog.md)
- Template workflow and command reference: [CLI](cli.md)
- API details and runtime option discovery: [API Reference](api_reference.md)
