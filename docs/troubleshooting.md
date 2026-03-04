# Troubleshooting

## `Unknown track alias`

Check supported aliases:

```python
from plotnado import Figure
Figure.available_track_aliases()
```

## Empty plots

- Confirm the region overlaps your data.
- Check chromosome naming consistency (`chr1` vs `1`).
- Validate DataFrame column names for your track type.

## TOML export/import errors

`to_toml()` requires `tomli-w`.

Install optional docs/dev extras:

```bash
pip install -e .[dev,docs]
```

## Matrix tracks fail to load

`cooler`-based tracks require cooler-compatible files and dependencies.

## Docs build issues

Reinstall docs extras and rebuild:

```bash
pip install -e .[docs]
mkdocs build --strict
```
