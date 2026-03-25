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

import typer

# Initialize the main CLI app
app = typer.Typer(
    help="PlotNado - Heuristic templates for genomic track visualization",
    no_args_is_help=True,
)

# Import commands - this registers them with the app
from . import init, plot, validate  # noqa: F401


def main():
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
