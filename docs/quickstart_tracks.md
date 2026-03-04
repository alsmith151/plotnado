# Build Tracks Fast

Use aliases with `Figure.add_track()` or chain helper methods.

## Alias workflow

```python
from plotnado import Figure

fig = Figure()
fig.add_track("scalebar")
fig.add_track("axis")
fig.add_track("genes", genome="hg38")
fig.add_track("bigwig", data="signal.bw", title="ChIP", color="#1f77b4", alpha=0.7)
```

## Shorthand kwargs

You can pass many aesthetics and label options directly:

- Aesthetics fields (for example `color`, `alpha`, `style`) are routed to `aesthetics`.
- Label fields (for example `plot_title`, `title_location`) are routed to `label`.

## Discover aliases and options

```python
from plotnado import Figure

Figure.available_track_aliases()
Figure.track_options("bigwig")
Figure.track_options_markdown("genes")
```

Or from CLI:

```bash
plotnado track-options
plotnado track-options bigwig
plotnado track-options --all --output-format json
```

Runnable example: `python examples/quickstart/02_aliases_and_options.py`

Example output:

![Alias and options quickstart](images/examples/quickstart_aliases_and_options.png)
