"""
Plotnado command-line interface.
"""

import pathlib
from enum import Enum
from typing import Optional

import typer
from typing_extensions import Annotated

import plotnado as pn


class OutputFormat(str, Enum):
    """Supported output formats."""

    png = "png"
    svg = "svg"
    pdf = "pdf"
    jpg = "jpg"
    jpeg = "jpeg"
    tiff = "tiff"


app = typer.Typer(help="Plotnado - Simple genomic track visualization")


@app.command()
def plot(
    coordinates: Annotated[
        str,
        typer.Argument(help="Coordinates to plot in format: CHR:START-END"),
    ],
    output: Annotated[
        Optional[pathlib.Path],
        typer.Option("--output", "-o", help="Output file path"),
    ] = None,
    output_format: Annotated[
        OutputFormat,
        typer.Option("--format", "-f", help="Output format"),
    ] = OutputFormat.png,
    width: Annotated[
        float,
        typer.Option("--width", "-w", help="Figure width in inches"),
    ] = 12.0,
    dpi: Annotated[
        int,
        typer.Option("--dpi", help="Resolution in dots per inch"),
    ] = 300,
):
    """
    Create a simple genome browser plot.

    Example:
        plotnado chr1:1000000-2000000 -o output.png
    """
    # Create a basic figure with scale and genes
    figure = pn.Figure(width=width)
    figure.add_track("scalebar")
    figure.add_track("genes", genome="hg38")

    # Determine output path
    if output is None:
        output = pathlib.Path(f"{coordinates.replace(':', '_')}.{output_format.value}")

    # Save the plot
    figure.save(output, coordinates, dpi=dpi)
    typer.echo(f"Saved plot to {output}")


def main():
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
