# Quick Start

This guide gets you from installation to a first saved figure quickly.

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

## 3) Next steps

- Learn all track-add patterns: [Track Construction](quickstart_tracks.md)
- Review practical defaults: [Best Practices](best_practices.md)
- Explore track families: [Track Catalog](track_catalog.md)
