import typer
import pathlib
from typing import Optional, List, Literal, Union
from typing_extensions import Annotated
from enum import Enum


class OuputFormat(str, Enum):
    png = "png"
    svg = "svg"
    pdf = "pdf"
    jpg = "jpg"
    jpeg = "jpeg"
    tiff = "tiff"
    bmp = "bmp"


def main():
    typer.run(plot)


def plot(
    template: Annotated[
        pathlib.Path,
        typer.Argument(
            exists=True,
            readable=True,
            resolve_path=True,
            help="Path to the template toml file.",
        ),
    ],
    coordinates: Annotated[
        Optional[str],
        typer.Option(
            "--coordinates",
            "-c",
            help="Coordinates to plot the template at. Provide in the format: CHR:START-END",
        ),
    ] = None,
    regions: Annotated[
        Optional[pathlib.Path],
        typer.Option("-r", "--regions", exists=True, readable=True, resolve_path=True),
    ] = None,
    name: Annotated[
        Optional[str],
        typer.Option(
            help="Name of the region to plot. Only used if regions is provided."
        ),
    ] = None,
    output_prefix: Optional[pathlib.Path] = pathlib.Path("plotnado_figures"),
    output_format: Annotated[
        OuputFormat,
        typer.Option("--format", "-f", help="Output format of the image."),
    ] = "png",
):
    """
    Use Plotnado to plot a template at either a set of coordinates or a set of regions.

    Args:
        template (pathlib.Path): Path to the template file.
        coordinates (str): Coordinates to plot the template at.
        regions (pathlib.Path): Path to the regions file.
        output (pathlib.Path): Path to the output image.
    """

    assert coordinates or regions, "Either coordinates or regions must be provided."

    import plotnado.api as pn

    figure = pn.Figure.from_toml(template)

    figures = dict()
    if coordinates:
        fig_name = f"{coordinates}" if not name else name
        figures[fig_name] = figure.plot(coordinates)
    elif regions:
        figures = figure.plot_regions(regions)

    for name, fig in figures.items():
        output_prefix.mkdir(exist_ok=True, parents=True)
        fmt = output_format.value
        fig.savefig(output_prefix / f"{name}.{fmt}", dpi=300)


if __name__ == "__main__":
    main()