# Installation

## Requirements

- Python 3.12+
- `uv` recommended, `pip` supported

## Install `uv`

If you do not already have `uv`, install it first:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

On macOS you can also use:

```bash
brew install uv
```

## Install PlotNado with `uv`

For a local virtual environment:

```bash
uv venv
source .venv/bin/activate
uv pip install plotnado
```

For a global CLI install:

```bash
uv tool install plotnado
```

## Install PlotNado with `pip`

If you prefer the standard Python packaging workflow:

```bash
python -m venv .venv
source .venv/bin/activate
pip install plotnado
```

## Install PlotNado in a Conda environment

If you already use Conda or Miniforge for scientific Python work:

```bash
conda create -n plotnado python=3.12
conda activate plotnado
pip install plotnado
```

For development in a Conda environment:

```bash
git clone https://github.com/alsmith151/plotnado
cd plotnado
conda create -n plotnado-dev python=3.12
conda activate plotnado-dev
pip install -e .[dev,docs]
pre-commit install
```

## Development install with `uv`

```bash
git clone https://github.com/alsmith151/plotnado
cd plotnado
uv venv
source .venv/bin/activate
uv sync --extra dev --extra docs
uv run pre-commit install
```

## Verify installation

```bash
uv run python -c "import plotnado; print(plotnado.__version__)"
uv run plotnado validate --help
```

## Build docs locally

```bash
uv run mkdocs build --strict
uv run mkdocs serve
```

## Run tests

```bash
uv run pytest tests/
```

## Development install with `pip`

If you prefer not to use `uv`, the traditional workflow still works:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .[dev,docs]
pre-commit install
```
