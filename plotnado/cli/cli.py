"""
Plotnado command-line interface.
"""

import json
import pathlib
from enum import Enum

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


def _markdown_cell(value: object) -> str:
    return str(value).replace("|", "\\|")


def _emit_options_table(track_alias: str, options: dict[str, dict], section: str | None) -> None:
    sections = [section] if section else ["track", "aesthetics", "label"]
    typer.echo(f"\n[{track_alias}]")
    for section_name in sections:
        section_options = options.get(section_name, {})
        typer.echo(f"  {section_name}:")
        if not section_options:
            typer.echo("    (none)")
            continue
        for field_name, meta in section_options.items():
            choices = meta.get("choices") or []
            choices_text = ",".join(str(choice) for choice in choices) if choices else "—"
            typer.echo(
                "    "
                f"{field_name}: type={meta['type']}, default={meta['default']}, "
                f"choices={choices_text}, required={meta['required']}"
            )


@app.command("track-options")
def track_options(
    track: Annotated[
        str | None,
        typer.Argument(
            help="Track alias to inspect (e.g. bigwig, genes, axis). Omit to list aliases."
        ),
    ] = None,
    all_tracks: Annotated[
        bool,
        typer.Option("--all", help="Show options for all track aliases"),
    ] = False,
    output_format: Annotated[
        str,
        typer.Option("--output-format", "-f", help="Output format: table, markdown, or json"),
    ] = "table",
    section: Annotated[
        str | None,
        typer.Option("--section", help="Optional section filter: track, aesthetics, or label"),
    ] = None,
):
    """Inspect available kwargs for each track alias.

    Examples:
        plotnado track-options
        plotnado track-options bigwig
        plotnado track-options bigwig --section aesthetics
        plotnado track-options bigwig -f markdown
        plotnado track-options --all -f json
    """
    aliases = pn.Figure.available_track_aliases()
    valid_sections = {"track", "aesthetics", "label"}
    valid_formats = {"table", "markdown", "json"}

    if output_format not in valid_formats:
        raise typer.BadParameter(
            f"--output-format must be one of: {', '.join(sorted(valid_formats))}"
        )
    if section is not None and section not in valid_sections:
        raise typer.BadParameter(f"--section must be one of: {', '.join(sorted(valid_sections))}")
    if track and all_tracks:
        raise typer.BadParameter("Use either a track alias argument or --all, not both")

    if track is None and not all_tracks:
        typer.echo("Available track aliases:")
        for alias, class_name in sorted(aliases.items()):
            typer.echo(f"  {alias:18} -> {class_name}")
        typer.echo("\nUse `plotnado track-options <alias>` for full option details.")
        return

    requested_aliases = [track] if track else sorted(aliases.keys())
    normalized_aliases = [alias.lower() for alias in requested_aliases]
    unknown_aliases = [alias for alias in normalized_aliases if alias not in aliases]
    if unknown_aliases:
        raise typer.BadParameter(
            f"Unknown alias(es): {', '.join(unknown_aliases)}. "
            f"Available: {', '.join(sorted(aliases.keys()))}"
        )

    if output_format == "json":
        payload = {
            alias: pn.Figure.track_options(alias)
            for alias in normalized_aliases
        }
        if section:
            payload = {
                alias: {section: data.get(section, {})}
                for alias, data in payload.items()
            }
        typer.echo(json.dumps(payload, indent=2))
        return

    if output_format == "markdown":
        for index, alias in enumerate(normalized_aliases):
            if section:
                options = pn.Figure.track_options(alias)
                typer.echo(f"## {alias}\n")
                typer.echo(f"### {section.title()} fields\n")
                typer.echo("| Name | Type | Default | Choices | Required | Description |")
                typer.echo("|---|---|---|---|---|---|")
                for field_name, meta in options.get(section, {}).items():
                    choices = meta.get("choices") or []
                    choices_text = _markdown_cell(
                        ", ".join(str(choice) for choice in choices) if choices else "—"
                    )
                    type_cell = _markdown_cell(meta["type"])
                    default_cell = _markdown_cell(meta["default"])
                    required_cell = _markdown_cell(meta["required"])
                    description_cell = _markdown_cell(meta.get("description") or "—")
                    typer.echo(
                        f"| {field_name} | {type_cell} | {default_cell} | {choices_text} | {required_cell} | {description_cell} |"
                    )
            else:
                typer.echo(pn.Figure.track_options_markdown(alias))
            if index != len(normalized_aliases) - 1:
                typer.echo("\n")
        return

    for alias in normalized_aliases:
        _emit_options_table(alias, pn.Figure.track_options(alias), section)


@app.command()
def plot(
    coordinates: Annotated[
        str,
        typer.Argument(help="Coordinates to plot in format: CHR:START-END"),
    ],
    output: Annotated[
        pathlib.Path | None,
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
