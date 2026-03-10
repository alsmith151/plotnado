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

## Docs build fails with Python version error

PlotNado requires Python 3.12+.

```bash
python --version
mkdocs build --strict
```

## TOML export/import errors

Install optional dependencies:

```bash
pip install -e .[docs]
```
