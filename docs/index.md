# PlotNado

![PlotNado logo](images/Logo.jpeg)

**Publication-ready genome browser figures in Python.**

PlotNado composes tracks — BigWig signals, gene annotations, peak intervals, interaction arcs — into a single figure with one fluent API call per track.

```bash
pip install plotnado
```

## First figure

```python
import numpy as np
import pandas as pd
from plotnado import GenomicFigure

bins = np.arange(1_000_000, 1_100_000, 1_000)
signal = pd.DataFrame({
    "chrom": "chr1", "start": bins,
    "end": bins + 1_000,
    "value": 5 + 2 * np.sin(np.linspace(0, 6, len(bins))),
})

fig = GenomicFigure()
fig.scalebar()
fig.axis()
fig.genes("hg38", title="Genes")
fig.bigwig(signal, title="ChIP signal", style="fill", color="#1f77b4")
fig.save("output.png", region="chr1:1,010,000-1,080,000")
```

![Quickstart figure](images/quickstart.png)

## Two workflows

**Python API** — build figures in code, ideal for notebooks and pipelines.

**CLI** — turn genomic files into an editable YAML template, then render it.

```bash
plotnado init *.bw peaks.narrowpeak --auto --output template.yaml
plotnado plot template.yaml --region chr1:1,000,000-1,100,000 --output out.png
```

## Where to go next

| Goal | Page |
| --- | --- |
| All methods and their parameters | [Reference](reference.md) |
| Styling — colors, labels, scales, themes | [Aesthetics](aesthetics.md) |
| CLI commands and YAML template format | [CLI](cli.md) |
