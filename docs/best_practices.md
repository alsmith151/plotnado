# Best Practices

Use this checklist when building figures for analysis reports or publications.

## Build figures deterministically

- Pin file paths and genome build explicitly (`genome="hg38"` rather than implicit defaults).
- Use `fig.to_toml("figure.toml")` to version figure definitions with your analysis code.
- Rebuild figures from TOML in CI or notebooks to avoid configuration drift.

## Keep track construction readable

- Prefer 4 to 8 tracks per panel for first-pass interpretation.
- Group structural tracks first (`scalebar`, `axis`, `genes`), then signals, then annotations.
- Use explicit titles on quantitative tracks so exported figures are self-describing.

## Keep colors semantically consistent

- Call `gf.autocolor()` once near figure setup.
- Reuse `color_group` across related tracks (for example BED + BigWig for one sample).
- Avoid manually overriding colors unless you need a specific publication palette.

```python
gf = GenomicFigure(theme="publication")
gf.autocolor()

gf.bed("A.bigBed", title="A peaks", color_group="A")
gf.bigwig("A.bw", title="A signal", color_group="A")
gf.bed("B.bigBed", title="B peaks", color_group="B")
gf.bigwig("B.bw", title="B signal", color_group="B")
```

## Use aliases for speed, models for strictness

- Use alias helpers (`fig.bigwig(...)`, `fig.bed(...)`) for exploratory work.
- Use explicit track classes when validating inputs/types in larger pipelines.
- Inspect options at runtime before custom styling:

```python
from plotnado import GenomicFigure

GenomicFigure.track_options("bigwig")
GenomicFigure.track_options_markdown("genes")
```

## Validate input assumptions early

- Confirm chromosome naming is consistent (`chr1` vs `1`).
- Check region overlap before styling tweaks.
- Keep coordinate systems and genome versions consistent across all inputs.

## Use themes intentionally

- Start with defaults for publication-ready output.
- Apply a project-wide theme once (instead of per-track overrides) to reduce style drift.
- Only disable theme (`theme=None`) when you need full manual control.

## Automate quality checks

- Run tests and docs build before release:

```bash
uv run pytest tests/ -v
uv run mkdocs build --strict
```

- Prefer example scripts in `examples/` as regression fixtures for plotting behavior.
