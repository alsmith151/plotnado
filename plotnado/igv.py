"""Parse IGV session XML files into plotnado Template objects."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

from .template import GuideSpec, GroupSpec, Template, TrackSpec
from .tracks.enums import TrackType


def _igv_color_to_hex(color_str: str) -> str | None:
    """Convert IGV RGB string '255,0,0' to hex '#ff0000'."""
    if not color_str:
        return None
    try:
        parts = color_str.split(",")
        if len(parts) != 3:
            return None
        r, g, b = (int(v.strip()) for v in parts)
        return f"#{r:02x}{g:02x}{b:02x}"
    except (ValueError, AttributeError):
        return None


def _infer_track_type(path: str) -> TrackType:
    """Infer TrackType from file extension."""
    ext = Path(path).suffix.lower()
    return {
        ".bw": TrackType.BIGWIG,
        ".bigwig": TrackType.BIGWIG,
        ".bed": TrackType.BED,
        ".bedgraph": TrackType.BEDGRAPH,
        ".bg": TrackType.BEDGRAPH,
        ".narrowpeak": TrackType.NARROWPEAK,
    }.get(ext, TrackType.UNKNOWN)


_ANNOTATION_CLAZZ = ("FeatureTrack", "SequenceTrack")
_GENE_ID_PATTERNS = ("ncbiRefSeq", "refseq", "gencode", "ensembl", "knownGenes")


def _is_annotation_track(track_elem: ET.Element) -> bool:
    """Return True if the track is a gene/sequence annotation (not a data track)."""
    clazz = track_elem.get("clazz", "")
    if any(c in clazz for c in _ANNOTATION_CLAZZ):
        return True
    track_id = track_elem.get("id", "").lower()
    name = track_elem.get("name", "").lower()
    return any(p.lower() in track_id or p.lower() in name for p in _GENE_ID_PATTERNS)


@dataclass
class IgvSession:
    """Parsed IGV session: a ready-to-use Template plus session metadata."""

    template: Template
    locus: str | None = None
    genome: str | None = None


def parse_igv_session(path: str | Path) -> IgvSession:
    """Parse an IGV session XML file into an :class:`IgvSession`.

    Args:
        path: Path to the IGV session ``.xml`` file.

    Returns:
        :class:`IgvSession` containing the compiled :class:`Template`, the
        session's saved locus, and the genome assembly name.

    Example::

        session = parse_igv_session("igv_session.xml")
        fig = GenomicFigure.from_template(session.template)
        fig.plot(session.locus)
    """
    tree = ET.parse(path)
    root = tree.getroot()

    genome = root.get("genome")
    locus = root.get("locus")

    tracks: list[TrackSpec] = []
    has_genes = False
    # group_name -> list of track titles (for GroupSpec construction)
    groups_seen: dict[str, list[str]] = {}

    for track_elem in root.iter("Track"):
        if _is_annotation_track(track_elem):
            has_genes = True
            continue

        if track_elem.get("visible", "true").lower() == "false":
            continue

        track_id = track_elem.get("id", "")
        name = track_elem.get("name", "")

        title = Path(name).stem if name else None
        color = _igv_color_to_hex(track_elem.get("color", ""))
        track_type = _infer_track_type(track_id)

        # Normalize IGV pixel height to plotnado units (~26px = 1.0)
        igv_height = float(track_elem.get("height", "26"))
        height = round(igv_height / 26.0, 2)

        # Map autoscaleGroup integer → named group string
        group_name: str | None = None
        autoscale_group = track_elem.get("autoscaleGroup")
        if autoscale_group:
            group_name = f"igv_group_{autoscale_group}"
            groups_seen.setdefault(group_name, [])
            if title:
                groups_seen[group_name].append(title)

        options: dict = {}
        auto_scale = track_elem.get("autoScale", "true").lower()
        if auto_scale == "false":
            data_range = track_elem.find("DataRange")
            if data_range is not None:
                min_val = data_range.get("minimum")
                max_val = data_range.get("maximum")
                if min_val is not None:
                    options["min_value"] = float(min_val)
                if max_val is not None:
                    options["max_value"] = float(max_val)

        tracks.append(TrackSpec(
            path=track_id or None,
            type=track_type,
            title=title,
            color=color,
            height=height,
            group=group_name,
            options=options,
        ))

    # Only create GroupSpecs for groups that appear on >1 track
    group_specs = [
        GroupSpec(name=gname, tracks=titles, autoscale=True, autocolor=False)
        for gname, titles in groups_seen.items()
        if len(titles) > 1
    ]

    template = Template(
        genome=genome,
        tracks=tracks,
        groups=group_specs,
        guides=GuideSpec(genes=has_genes, axis=True, scalebar=True),
    )

    return IgvSession(template=template, locus=locus, genome=genome)
