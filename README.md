# PlotNado

PlotNado is a Python package that provides a simple interface to create genome browser plots. It's a wrapper around the CoolBox library with a focus on simplicity and ease of use. It is designed to be used in Jupyter notebooks, but can also be used as a CLI.

## Installation

You can install PlotNado using pip:

```bash
pip install git+https://github.com/alsmith151/plotnado
```

## Usage

Here's a simple example of how to use PlotNado to create a genome browser plot:

```python
import plotnado as pn

figure = pn.Figure()
figure.add_track('scale')
figure.add_track('genes', 'hg38')

figure.plot('chr1:1000-2000')
```

This will create a plot of the region chr1:1000-2000 on the hg38 genome with a scale track and a genes track.

## Documentation

For more information on how to use PlotNado, check out the [documentation](https://plotnado.readthedocs.io/).

