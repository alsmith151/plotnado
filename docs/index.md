# PlotNado

![PlotNado logo](assets/plotnado-logo.svg)

PlotNado is a Python library for clean, publication-ready genome browser figures with a fast, chainable API.

## Install

```bash
pip install plotnado
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

## Read by task

- New environment setup: [Installation](installation.md)
- First successful figure: [Quick Start](quickstart.md)
- Ways to add tracks: [Track Construction](quickstart_tracks.md)
- Production guidance: [Best Practices](best_practices.md)
- Track families and support matrix: [Track Catalog](track_catalog.md)
- CLI and runtime option discovery: [CLI](cli.md) and [API Reference](api_reference.md)
