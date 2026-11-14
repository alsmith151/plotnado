"""
Plotnado command-line interface.

Primary workflow:
  1. plotnado init <track_files...> --output template.yaml
     Generates a template from track files using inference heuristics
  
  2. Edit template.yaml as needed (optional)
  
  3. plotnado plot template.yaml --region chr:start-end [--output out.png]
     Renders the template for the specified region
  
  4. plotnado validate template.yaml
     Validates and explains the template
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _package_version

import typer

# Initialize the main CLI app
app = typer.Typer(
    help="PlotNado - Heuristic templates for genomic track visualization",
    no_args_is_help=True,
)


def get_version() -> str:
    """Return the installed plotnado version."""
    try:
        return _package_version("plotnado")
    except PackageNotFoundError:
        return "unknown"


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"plotnado {get_version()}")
        raise typer.Exit()


@app.callback()
def _root(
    version: bool = typer.Option(
        False,
        "--version",
        "-V",
        help="Show the plotnado version and exit.",
        callback=_version_callback,
        is_eager=True,
    ),
) -> None:
    """PlotNado - Heuristic templates for genomic track visualization."""


@app.command("version")
def version_command() -> None:
    """Show the plotnado version."""
    typer.echo(f"plotnado {get_version()}")


# Import commands - this registers them with the app
from . import init, plot, validate  # noqa: F401


def main():
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
