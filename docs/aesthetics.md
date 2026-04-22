# Aesthetics

Every track kwarg is routed automatically — aesthetics fields go into the style model, label fields go into the label model. You do not need to nest anything manually:

```python
fig.bigwig(
    "signal.bw",
    title="H3K27ac",        # label
    title_color="black",    # label
    style="fill",           # aesthetics
    color="#1f77b4",        # aesthetics
    alpha=0.8,              # aesthetics
)
```

---

## BigWig styles

`style` controls how the signal is drawn. The four options:

```python
fig.bigwig(signal, title="fill",     style="fill",     color="#1f77b4")
fig.bigwig(signal, title="fragment", style="fragment", color="#d62728")
fig.bigwig(signal, title="scatter",  style="scatter",  color="#2ca02c", scatter_point_size=10)
fig.bigwig(signal, title="std",      style="std",      color="#9467bd")
```

![BigWig styles](images/aesthetics/bigwig_styles.png)

---

## Color and alpha

`color` accepts any matplotlib color string (hex, named, rgb tuple). `alpha` controls opacity 0–1.

```python
fig.bigwig(signal, title="color=#1f77b4  alpha=1.0", style="fill", color="#1f77b4", alpha=1.0)
fig.bigwig(signal, title="color=#d62728  alpha=0.8", style="fill", color="#d62728", alpha=0.8)
fig.bigwig(signal, title="color=#2ca02c  alpha=0.5", style="fill", color="#2ca02c", alpha=0.5)
fig.bigwig(signal, title="color=#ff7f0e  alpha=0.2", style="fill", color="#ff7f0e", alpha=0.2)
```

![Color and alpha](images/aesthetics/color_alpha.png)

---

## Automatic coloring

Call `autocolor()` once and PlotNado assigns palette colors to all tracks. Use `color_group` to share a color across tracks that represent the same biological sample.

```python
fig = GenomicFigure()
fig.autocolor()
fig.bigwig(signal_a, title="Sample A signal", color_group="A", style="fill")
fig.bed(peaks_a,    title="Sample A peaks",  color_group="A", display="expanded")
fig.bigwig(signal_b, title="Sample B signal", color_group="B", style="fill")
fig.bed(peaks_b,    title="Sample B peaks",  color_group="B", display="expanded")
```

![Autocolor and color_group](images/aesthetics/autocolor.png)

---

## Label placement

`title_location` positions the track title: `"left"` (default), `"right"`, or `"center"`. Set `label_on_track=True` to draw labels inside the track panel instead of in the margin.

```python
fig.bigwig(signal, title="title_location='left'",   title_location="left",   label_on_track=False)
fig.bigwig(signal, title="title_location='right'",  title_location="right",  label_on_track=False)
fig.bigwig(signal, title="title_location='center'", title_location="center", label_on_track=False)
fig.bigwig(signal, title="label_on_track=True",     title_location="left",   label_on_track=True)
```

![Label placement](images/aesthetics/label_placement.png)

---

## Scale annotations and label boxes

`plot_scale=True` adds a y-range annotation. `scale_location` places it `"left"` or `"right"`. `label_box_enabled=True` draws a translucent box behind label text.

```python
fig.bigwig(
    signal,
    title="label box + scale right",
    label_box_enabled=True,
    label_box_alpha=0.85,
    title_location="left",
    plot_scale=True,
    scale_location="right",
)
fig.bigwig(
    signal,
    title="no box, scale left",
    label_box_enabled=False,
    plot_scale=True,
    scale_location="left",
)
```

![Scale and label box](images/aesthetics/scale_box.png)

---

## BED and narrowPeak styling

BED tracks support `display="expanded"` for stacked intervals and `show_labels=True` for name text. narrowPeak adds score-based coloring via `color_by` and optional summit lines.

```python
fig.bed(
    bed_df,
    title="BED",
    display="expanded",
    show_labels=True,
    color="#4dac26",
)
fig.narrowpeak(
    narrowpeak_df,
    title="narrowPeak",
    color_by="signalValue",
    cmap="Oranges",
    show_labels=True,
    show_summit=True,
)
```

![BED and narrowPeak](images/aesthetics/bed_narrowpeak.png)

---

## Overlay, autoscale, and highlights

`autoscale(True)` + `autoscale_group` links y-axes across tracks. `overlay()` draws multiple signals in one panel. `highlight()` shades a region.

```python
fig = GenomicFigure()
fig.autoscale(True)
fig.highlight("chr1:1,032,000-1,046,000")
fig.highlight_style(color="#ffdd57", alpha=0.22)

fig.bigwig(signal_a, title="Sample A", autoscale_group="g1", style="fill", color="#1f77b4")
fig.bigwig(signal_b, title="Sample B", autoscale_group="g1", style="fill", color="#d62728")
fig.overlay([signal_c, signal_d], title="Overlay", colors=["#2ca02c", "#9467bd"], alpha=0.55)
```

![Overlay, autoscale, highlight](images/aesthetics/overlay_autoscale.png)

---

## Gene display modes

`display="collapsed"` fits all genes on one row. `display="expanded"` stacks overlapping genes.

```python
fig.genes("hg38", title="display='collapsed'", display="collapsed")
fig.genes("hg38", title="display='expanded'",  display="expanded")
```

![Gene display modes](images/aesthetics/gene_display.png)

---

## Links

Arc width and color can be driven by a `score` column using `color_by_score=True` and any matplotlib colormap.

```python
fig.links(links_df, title="Loops", color_by_score=True, cmap="viridis", alpha=0.8)
fig.hline(0.1, color="#666666", linewidth=0.8)
fig.vline(1_048_000, color="#222222", linewidth=1.0)
```

![Links, hline, vline](images/aesthetics/links.png)

---

## Themes

Pass a theme to `GenomicFigure` to apply a consistent style to all tracks. Built-in options: `"default"`, `"minimal"`, `"publication"`.

```python
fig = GenomicFigure(theme="publication")
```

For full manual control pass `theme=None`. Custom themes can be built with `plotnado.Theme`.

---

## Common label kwargs (all tracks)

| Kwarg | Default | Effect |
| --- | --- | --- |
| `title` | — | Track title text |
| `title_location` | `"left"` | `"left"` · `"right"` · `"center"` |
| `title_color` | `"#333333"` | Title text color |
| `title_size` | `10` | Title font size |
| `title_weight` | `"bold"` | `"bold"` · `"normal"` |
| `label_on_track` | `False` | Draw title inside the plot area |
| `label_box_enabled` | `True` | Background box behind labels |
| `label_box_alpha` | `0.9` | Label box opacity |
| `plot_scale` | `True` | Show y-range annotation |
| `scale_location` | `"right"` | `"left"` · `"right"` |
| `plot_title` | `True` | Toggle title visibility |
