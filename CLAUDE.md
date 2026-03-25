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

## Adding a New Track Type

1. Create `plotnado/tracks/mytrack.py` with a class extending `Track`
2. If it has a single primary color: `class MyAesthetics(BaseAesthetics)` (from `tracks/aesthetics.py`)
   If it renders multiple series: `class MyAesthetics(BaseMultiColorAesthetics)`
3. Register the track using the decorator:
   ```python
   from .registry import registry
   from .enums import TrackType

   @registry.register(TrackType.MYTRACK, aliases=["my_alias"])
   class MyTrack(Track): ...
   ```
4. Add `MYTRACK = "mytrack"` to `TrackType` in `tracks/enums.py`
5. Add a builder method to `GenomicFigure` in `figure_methods.py`:
   ```python
   def mytrack(self, data: Any, /, **kwargs) -> Self:
       return self.add_track(TrackType.MYTRACK, data=data, **kwargs)
   ```
6. Add a `MytrackKwargs` TypedDict entry and regenerate `_kwargs.py`:
   ```bash
   python scripts/generate_kwargs.py
   ```
7. Export from `plotnado/tracks/__init__.py` and `plotnado/__init__.py`
8. Write tests

## Method Map (render.py)

`TemplateCompiler` and `GenomicFigure` no longer maintain a separate method map.
Track lookup is centralized in `plotnado/tracks/registry.py`, and aliases such as
`bedgraph`, `annotation`, `unknown`, `bigwig_overlay`, and `scale` resolve there.

## Template Compilation Rules

- `TemplateCompiler.compile()` must **never mutate** the `Template` argument
- Resolved group indices go into `RenderPlan.resolved_group_indices`, not back into the template
- Group references are resolved case-insensitively against track `name` or `title` fields

## Common Pitfalls

- **Width override**: always use `width if width is not None else plan.width`, never `width or plan.width` (breaks when width=0.0)
- **TrackType source of truth**: use `plotnado.tracks.enums.TrackType` for both Python and template codepaths
- **Registry lookup**: if you need alias → class resolution, use `plotnado.tracks.registry.registry`, not a hard-coded method map
- **CLI shim**: `plotnado/cli/render.py` re-exports from `plotnado.render` — always import from `plotnado.render` in new code

## Testing

Run: `uv run pytest tests/`

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
uv venv
source .venv/bin/activate
uv pip install -e ".[dev]"
uv run pytest tests/
```

Entry point: `plotnado = "plotnado.cli.cli:main"` (defined in pyproject.toml)
