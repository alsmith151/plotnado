"""
CLI command: plotnado plot

Renders a template for specified genomic regions.
"""

from pathlib import Path
from typing import Optional
import json
import importlib.resources
import pandas as pd

import typer
from typing_extensions import Annotated
from rich.console import Console

import plotnado as pn
from plotnado.template import Template
from plotnado.tracks import GenomicRegion
from plotnado.cli.render import TemplateCompiler

from . import cli

console = Console()


def resolve_gene_region(gene_name: str, genome: Optional[str] = None) -> GenomicRegion:
    """
    Resolve a gene name to a genomic region.
    
    Args:
        gene_name: Gene symbol to look up (e.g., 'GNAQ')
        genome: Optional genome assembly (e.g., 'hg38', 'mm10')
    
    Returns:
        GenomicRegion corresponding to the gene
    
    Raises:
        ValueError: If gene not found or genome not available
    """
    if not genome:
        raise ValueError("Cannot resolve gene name without genome specification. Ensure template has genome defined.")
    
    # Load gene annotations from bundled data
    try:
        bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
        with open(bed_prefix / "genes.json") as handle:
            mapping = json.load(handle)
        
        if genome not in mapping:
            raise ValueError(f"Gene annotations not available for genome '{genome}'")
        
        gene_file = bed_prefix / mapping[genome]
        genes_df = pd.read_csv(gene_file, sep="\t", header=None)
    except Exception as e:
        raise ValueError(f"Failed to load gene annotations: {e}")
    
    # Parse BED format (chrom, start, end, name, ...)
    genes_df.columns = [
        "chrom",
        "start",
        "end",
        "name",
        *[f"field_{i}" for i in range(max(0, genes_df.shape[1] - 4))],
    ]
    
    # Match gene name (case-insensitive)
    match = genes_df.loc[genes_df["name"].astype(str).str.upper() == gene_name.upper()]
    
    if match.empty:
        raise ValueError(f"Gene '{gene_name}' not found in {genome} annotations")
    
    row = match.iloc[0]
    return GenomicRegion(
        chromosome=row["chrom"],
        start=int(row["start"]),
        end=int(row["end"]),
    )


@cli.app.command("plot")
def plot_command(
    template_file: Annotated[
        str,
        typer.Argument(help="Path to YAML template file"),
    ],
    region: Annotated[
        list[str],
        typer.Option(
            "--region",
            "-r",
            help="Genomic region to plot (chr:start-end or gene name). Repeat for multiple regions.",
        ),
    ],
    output: Annotated[
        Optional[str],
        typer.Option(
            "--output",
            "-o",
            help="Output image file path. Only valid with a single region.",
        ),
    ] = None,
    format: Annotated[
        str,
        typer.Option(
            "--format",
            "-f",
            help="Output format (png, pdf, svg, jpg)",
        ),
    ] = "png",
    width: Annotated[
        float,
        typer.Option(
            "--width",
            "-w",
            help="Figure width in inches",
        ),
    ] = 12.0,
    dpi: Annotated[
        int,
        typer.Option(
            "--dpi",
            help="Resolution (dots per inch)",
        ),
    ] = 600,
):
    """
    Render a template for one or more genomic regions.

    The plot command loads a template (created by 'plotnado init'),
    applies it to the specified region(s), and saves the resulting plot(s).

    Supports both genomic coordinates and gene names for the region parameter.

    Examples:
        plotnado plot template.yaml --region chr1:1000-2000
        plotnado plot template.yaml --region GNAQ
        plotnado plot template.yaml -r chr1:1M-2M -r chr2:5M-6M
        plotnado plot template.yaml -r chr1:1,000,000-2,000,000 -o plot.pdf
        plotnado plot template.yaml --region chr1:start-end --format svg --width 15
    """

    if not region:
        console.print("[red]Error: At least one --region is required[/red]")
        raise typer.Exit(code=1)

    if output and len(region) > 1:
        console.print("[red]Error: --output can only be used with a single --region[/red]")
        raise typer.Exit(code=1)

    # Load template
    try:
        template = Template.load(template_file)
        console.print(f"[green]✓ Loaded template:[/green] {template_file}")
    except FileNotFoundError:
        console.print(f"[red]Error: Template file not found:[/red] {template_file}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[red]Error loading template:[/red] {e}")
        raise typer.Exit(code=1)

    # Compile template once (region-independent)
    try:
        plan = TemplateCompiler.compile(template)
        console.print(f"[cyan]Compiled render plan:[/cyan] {len(plan.tracks)} tracks")
    except Exception as e:
        console.print(f"[red]Error compiling template:[/red] {e}")
        raise typer.Exit(code=1)

    for region_str in region:
        # Parse region - try gene name first if it doesn't look like a genomic coordinate
        gr = None

        if ":" not in region_str:
            try:
                gr = resolve_gene_region(region_str, genome=template.genome)
                console.print(f"[cyan]Resolved gene:[/cyan] {region_str} → {gr}")
            except ValueError as e:
                console.print(f"[yellow]Could not resolve as gene name:[/yellow] {e}")
                console.print("[yellow]Attempting to parse as genomic region...[/yellow]")

        if gr is None:
            try:
                gr = GenomicRegion.from_str(region_str)
                console.print(f"[cyan]Region:[/cyan] {gr}")
            except Exception as e:
                console.print(f"[red]Error parsing region '{region_str}':[/red] {e}")
                console.print("Expected format: chr:start-end (e.g., chr1:1000-2000) or gene name (e.g., GNAQ)")
                raise typer.Exit(code=1)

        # Create figure from render plan
        try:
            fig = pn.GenomicFigure(width=width or plan.width, track_height=plan.track_height)

            if plan.add_scalebar:
                fig.scalebar()
            if plan.add_axis:
                fig.axis()
            if plan.add_genes and plan.genome:
                fig.genes(plan.genome)

            for resolved_track in plan.tracks:
                method, data, kwargs = plan.get_track_by_method(resolved_track.index)
                track_method = getattr(fig, method)
                if data:
                    track_method(data, **kwargs)
                else:
                    track_method(**kwargs)
        except Exception as e:
            console.print(f"[red]Error creating figure for {gr}:[/red] {e}")
            raise typer.Exit(code=1)

        # Determine output path
        if output:
            out_path = Path(output)
        else:
            safe_region = str(gr).replace(":", "_").replace("-", "_").replace("(", "_").replace(")", "_")
            out_path = Path(f"{Path(template_file).stem}_{safe_region}.{format}")

        try:
            fig.save(out_path, region=str(gr), dpi=dpi)
            console.print(f"[green]✓ Saved plot:[/green] [bold]{out_path}[/bold]")
        except Exception as e:
            console.print(f"[red]Error saving figure:[/red] {e}")
            raise typer.Exit(code=1)
