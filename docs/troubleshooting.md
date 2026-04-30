# Troubleshooting

## `Unknown track alias`

```python
from plotnado import GenomicFigure

GenomicFigure.available_track_aliases()
```

## Empty plots

- Confirm the region overlaps your data.
- Check chromosome naming consistency (`chr1` vs `1`).
- Validate required columns for your track input.

## Overlay panel looks squashed or uses the wrong scale

Check these points in order:

- Put `autoscale_group` on the overlay track itself when the whole panel should match sibling signal tracks.
- Use `fig.autoscale(True)` if you want figure-wide shared limits derived from all signal-like tracks, including overlays.
- Remember that explicit overlay `min_value` / `max_value` intentionally overrides grouped or global autoscale on that edge.

```python
fig.overlay(
	[signal_a, signal_b],
	title="Overlay",
	autoscale_group="signal-panels",
)
```

If two series need fundamentally different scales to remain interpretable, split them into separate panels instead of forcing them into one overlay.

## Overlay label shows an unexpected range

The scale label mirrors the limits that were active when the overlay was rendered.

- If the label should match sibling panels, use `autoscale_group` on the overlay.
- If the label should stay pinned, set overlay `min_value` / `max_value` explicitly.

## Docs build fails with Python version error

PlotNado requires Python 3.12+.

```bash
uv run python --version
uv run mkdocs build --strict
```

## TOML export/import errors

Install optional dependencies:

```bash
uv sync --extra docs
```

## I cannot find a track option in the docs

Use runtime discovery first:

```python
from plotnado import GenomicFigure

GenomicFigure.track_options("overlay")
GenomicFigure.track_options_markdown("bigwig")
```

Then compare with [Aesthetics Reference](aesthetics_reference.md) and [Reference](reference.md).
