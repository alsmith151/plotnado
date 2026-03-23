"""
CLI command: plotnado init

Generates a template from track files using inference heuristics.
"""

from pathlib import Path
from typing import Optional

import typer
from typing_extensions import Annotated
from rich.console import Console
from rich.table import Table

from plotnado.template import Template, TrackSpec, GuideSpec, GroupSpec, TrackType
from plotnado.cli.inference import infer_track
from plotnado.cli.grouping import (
    PredefinedGroupingStrategies,
    apply_grouping_strategy,
    detect_and_apply_grouping,
)

from . import cli

console = Console()

# Track type ordering for generated templates: signal first, then peaks, then annotations
_TYPE_ORDER: dict[str, int] = {
    TrackType.BIGWIG.value: 0,
    TrackType.BEDGRAPH.value: 1,
    TrackType.NARROWPEAK.value: 2,
    TrackType.BED.value: 3,
    TrackType.ANNOTATION.value: 4,
    TrackType.LINKS.value: 5,
    TrackType.GENE.value: 6,
    TrackType.OVERLAY.value: 7,
    TrackType.UNKNOWN.value: 8,
}


def _sort_tracks(tracks: list[TrackSpec]) -> list[TrackSpec]:
    """Sort tracks: BigWig/Bedgraph → NarrowPeak → BED → Links → Unknown."""
    return sorted(tracks, key=lambda t: _TYPE_ORDER.get(str(t.type), 8))


def _show_preview_table(tracks: list[TrackSpec]) -> None:
    """Show a rich table preview of inferred tracks."""
    console.print(f"\n[bold cyan]Inferred {len(tracks)} track(s) — preview:[/bold cyan]\n")
    table = Table(show_header=True, header_style="bold")
    table.add_column("#", style="dim", width=3)
    table.add_column("Title")
    table.add_column("Type", width=12)
    table.add_column("Group")
    table.add_column("Color", width=9)
    table.add_column("File")

    for i, track in enumerate(tracks, 1):
        color_str = track.color or "—"
        group_str = track.group or "—"
        file_str = Path(track.path).name if track.path else "—"
        table.add_row(
            str(i),
            track.title or "—",
            str(track.type),
            group_str,
            color_str,
            file_str,
        )
    console.print(table)


@cli.app.command("init")
def init_command(
    tracks: Annotated[
        list[str],
        typer.Argument(help="Path or URL to track files (BigWig, BED, etc.)"),
    ],
    output: Annotated[
        str,
        typer.Option(
            "--output",
            "-o",
            help="Output YAML template file path",
        ),
    ] = "template.yaml",
    genome: Annotated[
        Optional[str],
        typer.Option(
            "--genome",
            "-g",
            help="Default genome (e.g., hg38, mm10)",
        ),
    ] = None,
    group_by: Annotated[
        Optional[str],
        typer.Option(
            "--group-by",
            help=(
                "Grouping strategy: predefined name (sample, antibody) or regex pattern. "
                "Examples: --group-by sample, --group-by '([^_]+)_rep[0-9]'"
            ),
        ),
    ] = None,
    auto: Annotated[
        bool,
        typer.Option(
            "--auto",
            help="Generate template automatically without prompting",
        ),
    ] = False,
    no_genes: Annotated[
        bool,
        typer.Option(
            "--no-genes",
            help="Do not include gene track by default",
        ),
    ] = False,
):
    """
    Generate a template from track files using inference heuristics.

    The init command analyzes your track files, infers track types and grouping,
    then generates an editable YAML template for rendering plots.

    Supports flexible grouping strategies:
    - Predefined: sample, antibody (for seqnado SAMPLE_ANTIBODY.bw patterns)
    - Regex: custom patterns like '([^_]+)_rep[0-9]' to group by filename prefix

    Examples:
        plotnado init sample1.bw sample2.bw peaks.bed
        plotnado init sample1_H3K27ac.bw sample1_H3K4me3.bw sample2_H3K27ac.bw \\
          --group-by sample
        plotnado init control_r1.bw control_r2.bw treat_r1.bw treat_r2.bw \\
          --group-by '([^_]+)_r[0-9]'
        plotnado init --auto *.bw
    """

    if not tracks:
        console.print("[red]Error: At least one track file is required[/red]")
        raise typer.Exit(code=1)

    console.print(f"\n[bold]PlotNado Template Generator[/bold]")
    console.print(f"Analyzing {len(tracks)} track file(s)...\n")

    # Run inference on all tracks
    template = Template()
    template.genome = genome
    template.guides = GuideSpec(
        genes=not no_genes,
        axis=True,
        scalebar=True,
    )

    inferences = []
    for track_path in tracks:
        result = infer_track(track_path)
        inferences.append((track_path, result))

        track_spec = TrackSpec(
            path=track_path,
            type=result.track_type,
            title=result.title,
            color=result.suggested_color,
        )
        template.tracks.append(track_spec)

    # Detect seqnado pattern
    seqnado_results = [infer.is_seqnado for _, infer in inferences]
    all_seqnado = all(seqnado_results)

    # Handle grouping
    grouping_result = None

    if group_by:
        # User provided explicit grouping strategy
        try:
            strategy = PredefinedGroupingStrategies.parse_group_by(group_by)
            if strategy:
                grouping_result = apply_grouping_strategy(tracks, strategy)
                if not grouping_result:
                    console.print(
                        f"[yellow]⚠ Grouping strategy '{group_by}' "
                        f"did not match any tracks[/yellow]"
                    )
            else:
                console.print(f"[red]Error: Unknown grouping strategy '{group_by}'[/red]")
                raise typer.Exit(code=1)
        except ValueError as e:
            console.print(f"[red]Error: {e}[/red]")
            raise typer.Exit(code=1)

    elif not auto:
        # Interactive mode
        if all_seqnado:
            # Seqnado files detected - offer seqnado strategies
            console.print("[bold cyan]Seqnado Pipeline Detected[/bold cyan]")

            samples = set()
            antibodies = set()
            for _, infer in inferences:
                if infer.is_seqnado:
                    samples.add(infer.seqnado_sample)
                    antibodies.add(infer.seqnado_antibody)

            console.print(f"  Samples: {', '.join(sorted(samples))}")
            console.print(f"  Antibodies: {', '.join(sorted(antibodies))}")

            console.print("\n[bold]How would you like to group tracks?[/bold]")
            console.print("  1. by sample (each antibody for same sample shares scaling)")
            console.print("  2. by antibody (each sample for same antibody shares scaling)")
            console.print("  3. no grouping")
            console.print("  4. custom regex pattern")

            choice = typer.prompt("Select option", type=int, default=1)

            if choice == 1:
                strategy = PredefinedGroupingStrategies.get("sample")
                grouping_result = apply_grouping_strategy(tracks, strategy)
            elif choice == 2:
                strategy = PredefinedGroupingStrategies.get("antibody")
                grouping_result = apply_grouping_strategy(tracks, strategy)
            elif choice == 4:
                pattern = typer.prompt(
                    "Enter regex pattern (e.g., '([^_]+)_rep[0-9]')"
                )
                try:
                    strategy = PredefinedGroupingStrategies.parse_group_by(pattern)
                    grouping_result = apply_grouping_strategy(tracks, strategy)
                except ValueError as e:
                    console.print(f"[red]Error: {e}[/red]")
                    raise typer.Exit(code=1)

        else:
            # Non-seqnado interactive wizard
            # 1. Genome prompt
            if not genome:
                genome_input = typer.prompt(
                    "Genome assembly",
                    default="hg38",
                    prompt_suffix=" [hg38/mm10/none]: ",
                    show_default=False,
                )
                if genome_input and genome_input.lower() != "none":
                    template.genome = genome_input

            # 2. Gene track prompt (only if genome set)
            if template.genome and not no_genes:
                include_genes = typer.confirm(
                    "Include gene annotation track?", default=True
                )
                template.guides.genes = include_genes
            else:
                template.guides.genes = False

            # 3. Grouping
            auto_result = detect_and_apply_grouping(tracks)
            if auto_result:
                console.print(f"\n[bold cyan]Detected grouping:[/bold cyan] {auto_result.explanation}")
                apply = typer.confirm("Apply this grouping?", default=True)
                if apply:
                    grouping_result = auto_result
            else:
                apply_manual = typer.confirm(
                    "\nGroup any tracks together for shared autoscaling?", default=False
                )
                if apply_manual:
                    # Show numbered list
                    console.print("\nTracks:")
                    for i, track in enumerate(template.tracks, 1):
                        console.print(f"  {i}. {track.title} ({track.type.value})")

                    # Collect groups interactively
                    manual_groups: dict[str, list[int]] = {}
                    while True:
                        indices_str = typer.prompt(
                            "Enter track numbers to group (e.g. 1,3) or leave empty to finish",
                            default="",
                        )
                        if not indices_str.strip():
                            break
                        try:
                            indices = [int(x.strip()) - 1 for x in indices_str.split(",")]
                            if any(i < 0 or i >= len(template.tracks) for i in indices):
                                console.print("[yellow]⚠ Invalid track numbers, try again[/yellow]")
                                continue
                            group_name = typer.prompt("Group name", default=f"group{len(manual_groups) + 1}")
                            manual_groups[group_name] = indices
                        except ValueError:
                            console.print("[yellow]⚠ Enter comma-separated numbers[/yellow]")

                    if manual_groups:
                        from plotnado.cli.grouping import GroupingResult
                        grouping_result = GroupingResult(
                            groups=manual_groups,
                            strategy_name="manual",
                            explanation=f"Manually grouped {len(manual_groups)} group(s)",
                        )

            # 4. BigWig style
            bw_tracks = [t for t in template.tracks if str(t.type) in (TrackType.BIGWIG.value, TrackType.BEDGRAPH.value)]
            if bw_tracks:
                style_input = typer.prompt(
                    "BigWig display style",
                    default="fill",
                    prompt_suffix=" [fill/line]: ",
                    show_default=False,
                )
                if style_input in ("fill", "line"):
                    for track in bw_tracks:
                        track.style = style_input

    else:
        # Auto mode: auto-detect the best grouping
        grouping_result = detect_and_apply_grouping(tracks)

    # Apply grouping result to template
    if grouping_result:
        console.print(f"\n[green]✓ {grouping_result.explanation}[/green]")

        for group_name, indices in grouping_result.groups.items():
            group_spec = GroupSpec(
                name=group_name,
                tracks=[template.tracks[i].title for i in indices],
                autoscale=True,
                autocolor=True,
            )
            template.groups.append(group_spec)

            # Set group on individual tracks
            for idx in indices:
                template.tracks[idx].group = group_name

        for group in template.groups:
            console.print(f"  {group.name}: {', '.join(str(t) for t in group.tracks)}")

    # Sort tracks: signal first, then peaks, then annotations
    template.tracks = _sort_tracks(template.tracks)

    # Show preview table
    _show_preview_table(template.tracks)

    # Save template with annotated header
    output_path = Path(output)
    args_str = " ".join(Path(t).name for t in tracks)
    template.save(output_path, header_args=args_str)

    console.print(f"[green]✓ Template saved to:[/green] [bold]{output_path}[/bold]\n")
    console.print("[bold]Next steps:[/bold]")
    console.print(f"  1. Review: cat {output_path}")
    console.print("  2. Edit as needed (optional)")
    console.print(f"  3. Plot: plotnado plot {output_path} --region chr1:1000-2000")
