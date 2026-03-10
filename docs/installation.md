# Installation

## Requirements

- Python 3.12+
- `pip`

## Install from PyPI

```bash
pip install plotnado
```

## Development install

```bash
git clone https://github.com/alsmith151/plotnado
cd plotnado
pip install -e .[dev,docs]
pre-commit install
```

## Verify installation

```bash
python -c "import plotnado; print(plotnado.__version__)"
plotnado track-options bigwig
```

## Build docs locally

```bash
mkdocs build --strict
mkdocs serve
```
