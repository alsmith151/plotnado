# CLAUDE.md — AI Coding Guidelines for plotnado

## What is plotnado?

plotnado is a lightweight Python package for creating genome browser-style plots from genomic data files (BigWig, BED, NarrowPeak, etc.). It provides two complementary interfaces:

1. **Python API** — fluent builder for programmatic use
2. **CLI + YAML templates** — file-driven workflow for non-programmers

## Architecture

```
plotnado/
├── figure.py          # GenomicFigure — the core plotting API
├── template.py        # Pydantic models for YAML templates (Template, TrackSpec, ...)
├── render.py          # TemplateCompiler → RenderPlan (template-to-figure bridge)
├── theme.py           # Theme / BuiltinTheme
├── tracks/            # Track implementations (BigWigTrack, BedTrack, ...)
│   ├── enums.py       # Internal enums incl. TrackType (figure-layer)
│   └── ...
└── cli/
    ├── cli.py         # Typer app entry point
    ├── init.py        # `plotnado init` — infer template from files
    ├── plot.py        # `plotnado plot` — render template for regions
    ├── validate.py    # `plotnado validate` — check template validity
    ├── inference.py   # Heuristics: infer track type/title from filename
    └── grouping.py    # Grouping strategies for init command
```

## Key Design Decisions

### Python API (GenomicFigure)

The fluent builder pattern allows chaining:

```python
fig = (
    GenomicFigure()
    .bigwig("signal.bw", title="H3K27ac")
    .narrowpeak("peaks.narrowpeak")
    .genes("hg38")
    .axis()
    .scalebar()
)
fig.save("out.png", region="chr1:1000000-2000000")
```

- Each method (`.bigwig()`, `.bed()`, etc.) appends a track and returns `self`
- `from_template(path)` builds a figure from a YAML template

### Template / CLI layer

YAML templates are human-readable, version-controllable, and editable:

```yaml
genome: hg38
tracks:
  - path: signal.bw
    type: bigwig
    title: H3K27ac
    group: sample1
guides:
  genes: true
  axis: true
  scalebar: true
```

`TemplateCompiler.compile(template)` → `RenderPlan` → `GenomicFigure` calls.

## TrackType vs TemplateTrackType

There are two separate enums — do not confuse them:

| Enum | Location | Purpose |
|---|---|---|
| `TrackType` | `plotnado/tracks/enums.py` | Internal figure enum; values match `GenomicFigure` method names |
| `TemplateTrackType` | `plotnado/template.py` | User-facing YAML vocabulary; values appear in template files |

`TemplateTrackType` has more values (e.g. `annotation`, `unknown`) that map to existing figure methods. The mapping is defined in `RenderPlan.get_track_by_method()` in `render.py`.

`TrackType = TemplateTrackType` alias exists in `template.py` for backward compatibility.

## Adding a New Track Type

1. Create `plotnado/tracks/mytrack.py` with a class extending `Track`
2. Add aesthetics class if needed; register field names
3. Add a method to `GenomicFigure` in `figure.py` (e.g. `.mytrack(data, **kwargs)`)
4. Add the alias in `GenomicFigure._alias_map()`
5. Add a `TemplateTrackType.MYTRACK = "mytrack"` value in `template.py`
6. Add the mapping in `RenderPlan.get_track_by_method()` in `render.py`
7. Export from `plotnado/tracks/__init__.py` and `plotnado/__init__.py`
8. Write tests

## Method Map (render.py)

`RenderPlan.get_track_by_method()` maps `TemplateTrackType` values to `GenomicFigure` method names:

| TemplateTrackType | GenomicFigure method | Notes |
|---|---|---|
| `bigwig` | `bigwig` | |
| `bedgraph` | `bigwig` | BigWigTrack handles bedgraph natively |
| `bed` | `bed` | |
| `narrowpeak` | `narrowpeak` | |
| `gene` | `genes` | |
| `links` | `links` | |
| `annotation` | `bed` | BED interval track with annotation semantics |
| `overlay` | `overlay` | |
| `unknown` | `bed` | Fallback |

## Template Compilation Rules

- `TemplateCompiler.compile()` must **never mutate** the `Template` argument
- Resolved group indices go into `RenderPlan.resolved_group_indices`, not back into the template
- Group references are resolved case-insensitively against track `name` or `title` fields

## Common Pitfalls

- **Width override**: always use `width if width is not None else plan.width`, never `width or plan.width` (breaks when width=0.0)
- **bedgraph method**: there is no `GenomicFigure.bedgraph()` — bedgraph files use `.bigwig()`
- **TemplateTrackType vs TrackType**: don't import the wrong one; `template.py` is for template context, `tracks/enums.py` is for figure context
- **CLI shim**: `plotnado/cli/render.py` re-exports from `plotnado.render` — always import from `plotnado.render` in new code

## Testing

Run: `pytest tests/`

| Test file | Coverage |
|---|---|
| `test_template.py` | Template round-trip, YAML serialization |
| `test_render.py` | TemplateCompiler, no-mutation guarantee, autocolor |
| `test_inference.py` | Track type/title inference heuristics |
| `test_grouping.py` | Grouping strategies |
| `test_cli.py` | CLI integration via `typer.testing.CliRunner` |

Guidelines:
- Use `tmp_path` pytest fixture for file-based tests
- Do not mock `GenomicFigure.plot()` or `.save()` in unit tests
- The no-mutation guarantee on `TemplateCompiler.compile()` must always have a regression test

## Dev Setup

```bash
source .venv/bin/activate
uv pip install -e ".[dev]"
pytest tests/
```

Entry point: `plotnado = "plotnado.cli.cli:main"` (defined in pyproject.toml)
