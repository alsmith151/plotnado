# CLI

PlotNado includes a lightweight CLI.

## Plot a region

```bash
plotnado plot chr1:1,000,000-1,100,000 -o output.png
```

## Discover track options

```bash
plotnado track-options
plotnado track-options bigwig
plotnado track-options bigwig --section aesthetics
plotnado track-options --all --output-format json
plotnado track-options genes --output-format markdown
```

The CLI option metadata is generated from the same Pydantic models used by the Python API.
