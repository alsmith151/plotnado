# Track Construction

This page shows the different ways to add tracks to a `GenomicFigure`.

## 1) Preferred chain-helper style

This is concise and works well for interactive analysis.

```python
from plotnado import GenomicFigure

gf = GenomicFigure(theme="publication")
gf.autocolor()
gf.scalebar()
gf.genes("hg38", display="expanded", minimum_gene_length=1e5)

gf.bed(
    "https://.../THP1H3K4me3_bigBed.bigBed",
    title="THP1 H3K4me3 peaks",
    title_color="black",
    color_group="THP1 H3K4me3",
    draw_edges=False,
)

gf.bigwig(
    "https://.../THP1H3K4me1_bigWig.bigWig",
    title="THP1 H3K4me1",
    title_color="black",
    color_group="THP1 H3K4me1",
    style="std",
)

gf.axis()
gf.plot_gene("GNAQ")
```

## 2) Generic alias entrypoint (`add_track`)

Best when alias names come from config/TOML/CLI inputs.

```python
gf.add_track("scalebar")
gf.add_track("genes", genome="hg38")
gf.add_track("bigwig", data="signal.bw", title="ChIP", style="std")
```

## 3) Explicit track objects

Best when you want strict object construction and reuse track instances.

```python
from plotnado import BigWigTrack

gf.add_track(
    BigWigTrack(
        data="signal.bw",
        title="ChIP",
        style="std",
    )
)
```

## 4) Mixed approach (recommended in real projects)

- Use helper methods for most tracks.
- Use explicit classes only for custom/advanced tracks.
- Keep structural tracks first (`scalebar`, `axis`, `genes`) for readability.

## Automatic kwarg routing in `GenomicFigure` methods

When you use helper methods like `gf.bigwig(...)`, PlotNado can automatically route kwargs into nested models:

- Track kwargs stay at track level (for example `data`, `autoscale_group`).
- Aesthetics kwargs are auto-packed into `aesthetics` (for example `style`, `color`, `alpha`).
- Label kwargs are auto-packed into `label` (for example `title`, `title_color`, `title_location`).

Shorthand (recommended):

```python
gf.bigwig(
    "signal.bw",
    title="H3K4me3",
    title_color="black",
    style="std",
    color="#1f77b4",
    alpha=0.8,
)
```

Equivalent explicit form:

```python
gf.bigwig(
    "signal.bw",
    aesthetics={"style": "std", "color": "#1f77b4", "alpha": 0.8},
    label={"title": "H3K4me3", "title_color": "black"},
)
```

## Color management (`autocolor`, `color_group`)

Use `autocolor()` to assign colors automatically from the active theme.

```python
gf = GenomicFigure(theme="publication")
gf.autocolor()
gf.bigwig("sample_a.bw", title="Sample A")
gf.bigwig("sample_b.bw", title="Sample B")
```

Use `color_group` when multiple tracks should share the same assigned color.

```python
gf = GenomicFigure(theme="publication")
gf.autocolor()

gf.bed("peaks_a.bigBed", title="A peaks", color_group="A")
gf.bigwig("signal_a.bw", title="A signal", color_group="A")

gf.bed("peaks_b.bigBed", title="B peaks", color_group="B")
gf.bigwig("signal_b.bw", title="B signal", color_group="B")
```

Practical rule:

- Same biological sample/condition: same `color_group`.
- Different sample/condition: different `color_group`.

## Option discovery

```python
from plotnado import GenomicFigure

GenomicFigure.available_track_aliases()
GenomicFigure.track_options("bigwig")
GenomicFigure.track_options_markdown("genes")
```

See [Track Aliases](track_aliases.md) for alias notes, [Example Coverage](example_coverage.md) for runnable script mapping, and [Best Practices](best_practices.md) for production guidance.
