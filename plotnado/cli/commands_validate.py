"""
CLI command: plotnado validate

Validates and explains a template.
"""

from pathlib import Path

import typer
from typing_extensions import Annotated
from rich.console import Console
from rich.table import Table

from plotnado.template import Template
from plotnado.cli.render import TemplateCompiler

from . import cli

console = Console()


@cli.app.command("validate")
def validate_command(
    template_file: Annotated[
        str,
        typer.Argument(help="Path to YAML template file"),
    ],
):
    """
    Validate a template and show its configuration.

    Checks for missing files, group reference errors, and other issues
    before you attempt to plot.

    Examples:
        plotnado validate template.yaml
    """

    # Load template
    try:
        template = Template.load(template_file)
        console.print(f"[green]✓ Template loaded:[/green] {template_file}\n")
    except FileNotFoundError:
        console.print(f"[red]Error: Template file not found:[/red] {template_file}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[red]Error loading template:[/red] {e}")
        raise typer.Exit(code=1)

    # Show metadata
    console.print("[bold cyan]Metadata[/bold cyan]")
    if template.genome:
        console.print(f"  genome: {template.genome}")
    console.print(f"  width: {template.width} inches")
    console.print(f"  track_height: {template.track_height}")

    # Show guides
    console.print("\n[bold cyan]Guides[/bold cyan]")
    console.print(f"  genes: {template.guides.genes}")
    console.print(f"  axis: {template.guides.axis}")
    console.print(f"  scalebar: {template.guides.scalebar}")

    # Show tracks
    console.print(f"\n[bold cyan]Tracks ({len(template.tracks)})[/bold cyan]")

    table = Table(show_header=True, header_style="bold")
    table.add_column("Index", style="dim")
    table.add_column("Title")
    table.add_column("Type")
    table.add_column("Group")
    table.add_column("File")

    missing_files: list[str] = []

    for i, track in enumerate(template.tracks, 1):
        group_str = track.group or "—"
        file_str = Path(track.path).name if track.path else "—"

        # Check file existence for local paths
        if track.path and not track.path.startswith(("http://", "https://", "s3://", "ftp://")):
            if not Path(track.path).exists():
                missing_files.append(track.path)
                file_str = f"[red]{file_str} ✗[/red]"

        table.add_row(
            str(i),
            track.title or "—",
            str(track.type),
            group_str,
            file_str,
        )

    console.print(table)

    # Show groups
    if template.groups:
        console.print(f"\n[bold cyan]Groups ({len(template.groups)})[/bold cyan]")
        for group in template.groups:
            console.print(f"  {group.name}")
            if group.tracks:
                console.print(f"    tracks: {group.tracks}")
            console.print(f"    autoscale: {group.autoscale}")
            console.print(f"    autocolor: {group.autocolor}")

    # Report missing files
    errors = False
    if missing_files:
        errors = True
        console.print(f"\n[red]✗ {len(missing_files)} file(s) not found:[/red]")
        for f in missing_files:
            console.print(f"    {f}")

    # Dry-compile to catch group reference errors
    try:
        TemplateCompiler.compile(template)
        console.print("\n[green]✓ Group references resolved successfully[/green]")
    except ValueError as e:
        errors = True
        console.print(f"\n[red]✗ Group reference error:[/red] {e}")
        console.print("[dim]Tip: Check that track titles in your groups section match track titles exactly (case-insensitive).[/dim]")

    if errors:
        console.print("\n[red]Validation failed — fix the issues above before plotting.[/red]")
        raise typer.Exit(code=1)

    console.print("\n[green]✓ Validation complete[/green]")
    console.print("Next: plotnado plot <template> --region chr:start-end")
