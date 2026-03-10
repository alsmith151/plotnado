# Quick Start

This guide gets you from install to a first figure quickly.

## 1) Create a minimal figure

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

`GenomicFigure()` applies publication-quality styling by default. Pass `theme=None` for unthemed rendering.

## 2) Run a ready-made script

```bash
python examples/quickstart/01_first_plot.py
```

Expected output:

![Quickstart first plot](images/examples/quickstart_first_plot.png)

## 3) Next steps

- Add tracks quickly with aliases: [Build Tracks Fast](quickstart_tracks.md)
- Explore track families: [Track Catalog](track_catalog.md)
- See practical patterns: [Recipes](recipes.md)
