# Reference

## GenomicFigure

```python
from plotnado import GenomicFigure

fig = GenomicFigure(
    width=12,           # figure width in inches
    track_height=1.25,  # height multiplier per track
    theme="publication" # "default" | "minimal" | "publication" | None
)
```

Every track method returns `self`, so calls can be chained or written sequentially.

---

## Track methods

### Signal

#### `bigwig(data, /, **kwargs)`

Continuous signal. `data` is a BigWig file path or a bedgraph-like DataFrame with columns `chrom`, `start`, `end`, `value`.

Key kwargs:

| Kwarg | Default | Description |
| --- | --- | --- |
| `style` | `"std"` | `"fill"` · `"fragment"` · `"scatter"` · `"std"` |
| `color` | `"steelblue"` | Fill/line color |
| `alpha` | `1.0` | Opacity 0–1 |
| `linewidth` | `1.0` | Stroke width |
| `scatter_point_size` | `1.0` | Marker size (scatter only) |
| `smoothing_window` | — | Bin width for rolling average |
| `smoothing_method` | `"mean"` | `"mean"` or `"median"` |
| `min_value` / `max_value` | — | Fixed y-axis limits |

#### `overlay(tracks, /, **kwargs)`

Multiple signals in one panel. `tracks` is a list of data sources (same formats as `bigwig`).

Key kwargs: `colors` (list), `alpha`, `linewidth`.

#### `bigwig_diff(data_a, data_b, /, **kwargs)`

Difference track between two BigWig sources. `method` selects `"subtraction"`, `"ratio"`, or `"log2ratio"`.

#### `bigwig_collection(files, /, **kwargs)`

Multi-BigWig panel rendered as a heatmap or stacked lines. Pass a list of file paths.

---

### Intervals

#### `bed(data, /, **kwargs)`

Rectangle intervals. `data` is a BED/BigBed file path or DataFrame with `chrom`, `start`, `end`, optional `name`.

Key kwargs:

| Kwarg | Default | Description |
| --- | --- | --- |
| `display` | `"collapsed"` | `"collapsed"` stacks on one row; `"expanded"` stacks overlaps |
| `color` | `"steelblue"` | Fill color |
| `alpha` | `1.0` | Opacity |
| `show_labels` | `False` | Show interval name labels |
| `label_field` | `"name"` | Column used for labels |
| `draw_edges` | `True` | Outline rectangles |
| `interval_height` | `0.45` | Rectangle height (normalized) |
| `max_rows` | `5` | Maximum stacked rows |

#### `narrowpeak(data, /, **kwargs)`

ENCODE narrowPeak format with summit markers and score-based coloring.

Key kwargs (in addition to `bed` kwargs):

| Kwarg | Default | Description |
| --- | --- | --- |
| `color_by` | — | Column name to map color from (e.g. `"signalValue"`) |
| `cmap` | `"Oranges"` | Colormap name when `color_by` is set |
| `show_summit` | `False` | Draw a vertical summit line |

#### `links(data, /, **kwargs)`

Interaction arcs (BEDPE format). `data` has columns `chrom1`, `start1`, `end1`, `chrom2`, `start2`, `end2`, optional `score`.

Key kwargs: `color_by_score` (bool), `cmap`, `alpha`, `linewidth`.

#### `highlights(data, /, **kwargs)`

Shaded region overlays from a BED file. Use `highlight(region_str)` for inline string regions:

```python
fig.highlight("chr1:1,032,000-1,046,000")
fig.highlight_style(color="#ffdd57", alpha=0.22)
```

#### `hline(y_value, /, **kwargs)` · `vline(x_position, /, **kwargs)`

Horizontal or vertical reference lines. `color` and `linewidth` are the main kwargs.

---

### Annotation

#### `genes(genome="hg38", /, **kwargs)`

Gene models from bundled annotations (`"hg38"`, `"mm10"`) or a custom BED12/GTF file path via `data=`.

Key kwargs:

| Kwarg | Default | Description |
| --- | --- | --- |
| `display` | `"collapsed"` | `"collapsed"` · `"expanded"` |
| `minimum_gene_length` | `0` | Hide genes shorter than this (bp) |
| `max_number_of_rows` | `4` | Row limit in expanded mode |
| `show_labels` | `True` | Gene name text |
| `gene_label_font_size` | `10` | Label font size |
| `color` | `"steelblue"` | Exon fill color |

#### `scalebar(**kwargs)`

Genomic distance bar. Set `plot_scale=True` to show the numeric distance. Use `scale_location` (`"left"` · `"right"`) to position the text.

#### `axis(**kwargs)`

Coordinate tick labels. `num_ticks` controls tick density.

#### `spacer(height=0.5, **kwargs)`

Blank vertical gap between tracks.

---

## Figure-level controls

```python
fig.autoscale(True)                # link y-axes of tracks that share autoscale_group
fig.autocolor()                    # assign palette colors to tracks automatically
fig.autocolor("tab10")             # use a specific matplotlib palette
fig.highlight("chr1:1000000-1100000")     # shade a region on all tracks
fig.highlight_style(color="#ffdd57", alpha=0.22)
fig.defaults(genome="hg38")        # shortcut for autocolor + genes defaults
```

### Autoscale groups

Tracks sharing the same `autoscale_group` string are scaled to the same y-axis range after all data is loaded:

```python
fig.autoscale(True)
fig.bigwig("a.bw", title="A", autoscale_group="chip")
fig.bigwig("b.bw", title="B", autoscale_group="chip")
```

### Color groups

When using `autocolor()`, tracks with the same `color_group` string share one palette slot:

```python
fig.autocolor()
fig.bigwig("a.bw", title="A signal", color_group="sample_a")
fig.bed("a.bed",   title="A peaks",  color_group="sample_a")
fig.bigwig("b.bw", title="B signal", color_group="sample_b")
```

---

## Plotting and saving

```python
fig.plot("chr1:1,000,000-2,000,000")                                         # returns matplotlib Figure
fig.plot("chr1:1,000,000-2,000,000", extend=0.25)                            # extend viewport 25% each side
fig.plot_gene("GNAQ", extend=0.5)                                            # center on a gene
fig.plot_regions(["chr1:1,000,000-2,000,000", "chr2:5,000,000-6,000,000"], ncols=2)
fig.save("out.png", region="chr1:1,000,000-2,000,000", dpi=200)
fig.save("out.pdf", region="GNAQ")
```

---

## Editing after build

```python
fig["track title"]          # access track by title (case-insensitive)
fig[0]                      # access by index
fig["track title"].color = "coral"  # mutate in place

fig.update_track("title", color="purple", max_value=500)   # one track
fig.update_track(height=0.2)                               # all tracks
fig.update_track(track_type="bigwig", height=0.3)          # by type
fig.update_track(group="group_1", max_value=1000)          # by autoscale group
fig.update_track(where=lambda t: "CD" in (t.title or ""), color="coral")

fig.remove_track("title")
fig.remove_track(0)

fig.bigwig("extra.bw", title="extra", position="top")     # prepend
```

All editing methods return `self`.

---

## Load and save figure configs

```python
# From a YAML template (CLI or hand-written)
fig = GenomicFigure.from_template("template.yaml")

# From an IGV session XML
fig, locus = GenomicFigure.from_igv_session("session.xml")
fig.plot(locus)

# TOML round-trip
fig.to_toml("figure.toml")
fig2 = GenomicFigure.from_toml("figure.toml")
```

---

## Runtime option discovery

```python
GenomicFigure.available_track_aliases()           # all registered aliases
GenomicFigure.track_options("bigwig")             # dict of all kwargs + defaults
GenomicFigure.track_options_markdown("bigwig")    # print-friendly table
```
