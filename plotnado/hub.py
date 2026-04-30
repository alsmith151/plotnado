"""Parse UCSC-style track hubs into plotnado templates."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import urljoin, urlparse
from urllib.request import urlopen

from .cli.inference import infer_track
from .template import Template, TrackSpec
from .tracks.enums import TrackType


_REMOTE_SCHEMES = {"http", "https", "ftp"}
_HIDDEN_STATES = {"hide", "off", "0"}
_SUPPORTED_TRACK_TYPES = {
    TrackType.BIGWIG,
    TrackType.BEDGRAPH,
    TrackType.BED,
    TrackType.NARROWPEAK,
    TrackType.LINKS,
}
_TRACK_TYPE_MAP = {
    "bigwig": TrackType.BIGWIG,
    "bedgraph": TrackType.BEDGRAPH,
    "bigbed": TrackType.BED,
    "bed": TrackType.BED,
    "narrowpeak": TrackType.NARROWPEAK,
    "links": TrackType.LINKS,
    "biginteract": TrackType.LINKS,
    "interact": TrackType.LINKS,
}


@dataclass
class UcscHubSession:
    """Parsed UCSC track hub ready for template-based rendering."""

    template: Template
    genome: str
    source: str
    short_label: str | None = None
    long_label: str | None = None


def parse_ucsc_hub(
    source: str | Path,
    *,
    genome: str | None = None,
    include_hidden: bool = True,
) -> UcscHubSession:
    """Parse a UCSC hub.txt entrypoint into a :class:`Template`.

    Args:
        source: Hub entrypoint as a local path or remote URL.
        genome: Optional genome id to select when the hub serves multiple genomes.
        include_hidden: Preserve hidden tracks as zero-height tracks.

    Returns:
        Parsed hub session with a plotnado template.
    """

    source_str = str(source)
    hub_records = _parse_stanzas(_read_text(source_str))
    if not hub_records:
        raise ValueError(f"Hub source {source_str!r} did not contain any records")

    hub_record = hub_records[0]
    genomes_ref = hub_record.get("genomesFile")
    if not genomes_ref:
        raise ValueError(f"Hub source {source_str!r} is missing a genomesFile entry")

    genomes_source = _resolve_reference(source_str, genomes_ref)
    genome_records = _parse_stanzas(_read_text(genomes_source))
    if not genome_records:
        raise ValueError(f"Genome source {genomes_source!r} did not contain any records")

    genome_record = _select_genome_record(genome_records, genome)
    track_db_ref = genome_record.get("trackDb")
    if not track_db_ref:
        raise ValueError(
            f"Genome {genome_record.get('genome', '<unknown>')!r} is missing a trackDb entry"
        )

    track_db_source = _resolve_reference(genomes_source, track_db_ref)
    track_records = [record for record in _parse_stanzas(_read_text(track_db_source)) if "track" in record]
    track_specs = _build_track_specs(
        track_records,
        track_db_source=track_db_source,
        include_hidden=include_hidden,
    )

    template = Template(
        genome=genome_record["genome"],
        tracks=track_specs,
    )
    return UcscHubSession(
        template=template,
        genome=genome_record["genome"],
        source=source_str,
        short_label=hub_record.get("shortLabel"),
        long_label=hub_record.get("longLabel"),
    )


def _build_track_specs(
    records: Iterable[dict[str, str]],
    *,
    track_db_source: str,
    include_hidden: bool,
) -> list[TrackSpec]:
    record_list = list(records)
    records_by_id = {record["track"]: record for record in record_list}
    track_specs: list[TrackSpec] = []

    for record in record_list:
        data_ref = record.get("bigDataUrl")
        if not data_ref:
            continue

        resolved_path = _resolve_reference(track_db_source, data_ref)
        track_type = _infer_track_type(record, resolved_path)
        if track_type is None:
            continue

        hidden = _is_hidden(record, records_by_id)
        if hidden and not include_hidden:
            continue

        track_specs.append(
            TrackSpec(
                path=resolved_path,
                type=track_type,
                title=record.get("longLabel") or record.get("shortLabel") or record["track"],
                name=record["track"],
                color=_ucsc_color_to_hex(record.get("color")),
                height=0.0 if hidden else 1.0,
            )
        )

    return track_specs


def _infer_track_type(record: dict[str, str], resolved_path: str) -> TrackType | None:
    type_value = record.get("type", "").split()
    if type_value:
        mapped = _TRACK_TYPE_MAP.get(type_value[0].lower())
        if mapped in _SUPPORTED_TRACK_TYPES:
            return mapped

    inferred = infer_track(resolved_path).track_type
    if inferred in _SUPPORTED_TRACK_TYPES:
        return inferred
    return None


def _is_hidden(
    record: dict[str, str],
    records_by_id: dict[str, dict[str, str]],
    seen: set[str] | None = None,
) -> bool:
    track_id = record.get("track")
    if track_id is None:
        return False

    if seen is None:
        seen = set()
    if track_id in seen:
        return False
    seen.add(track_id)

    visibility = record.get("visibility", "").split()
    if visibility and visibility[0].lower() in _HIDDEN_STATES:
        return True

    parent_id, parent_state = _split_reference_state(record.get("parent"))
    if parent_state in _HIDDEN_STATES:
        return True
    if parent_id and parent_id in records_by_id and _is_hidden(records_by_id[parent_id], records_by_id, seen):
        return True

    super_id, super_state = _split_reference_state(record.get("superTrack"))
    if super_state in _HIDDEN_STATES:
        return True
    if super_id and super_id in records_by_id and _is_hidden(records_by_id[super_id], records_by_id, seen):
        return True

    return False


def _split_reference_state(value: str | None) -> tuple[str | None, str | None]:
    if not value:
        return None, None
    parts = value.split()
    track_id = parts[0] if parts else None
    state = parts[1].lower() if len(parts) > 1 else None
    return track_id, state


def _select_genome_record(records: list[dict[str, str]], genome: str | None) -> dict[str, str]:
    available = [record.get("genome") for record in records if record.get("genome")]
    if genome is None:
        if len(records) == 1 and records[0].get("genome"):
            return records[0]
        raise ValueError(f"Hub serves multiple genomes; choose one of {available}")

    for record in records:
        if record.get("genome") == genome:
            return record

    raise ValueError(f"Genome {genome!r} not found in hub; available genomes: {available}")


def _parse_stanzas(text: str) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] = {}
    last_key: str | None = None
    pending_key: str | None = None

    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            if pending_key is not None:
                current[pending_key] = ""
                last_key = pending_key
                pending_key = None
            if current:
                records.append(current)
                current = {}
                last_key = None
            continue

        if pending_key is not None:
            current[pending_key] = stripped
            last_key = pending_key
            pending_key = None
            continue

        parts = stripped.split(None, 1)
        if len(parts) == 1:
            pending_key = parts[0]
            continue

        key, value = parts[0], parts[1].strip()
        current[key] = value
        last_key = key

    if pending_key is not None:
        current[pending_key] = ""
    if current:
        records.append(current)

    return records


def _read_text(source: str) -> str:
    parsed = urlparse(source)
    if parsed.scheme in _REMOTE_SCHEMES:
        with urlopen(source) as response:  # nosec: trusted user-supplied data source
            return response.read().decode("utf-8")
    if parsed.scheme == "file":
        return Path(parsed.path).read_text(encoding="utf-8")
    return Path(source).read_text(encoding="utf-8")


def _resolve_reference(base: str, reference: str) -> str:
    parsed_reference = urlparse(reference)
    if parsed_reference.scheme in _REMOTE_SCHEMES or parsed_reference.scheme == "file":
        return reference

    parsed_base = urlparse(base)
    if parsed_base.scheme in _REMOTE_SCHEMES:
        return urljoin(base, reference)
    if parsed_base.scheme == "file":
        base_path = Path(parsed_base.path)
    else:
        base_path = Path(base)
    return str((base_path.parent / reference).resolve())


def _ucsc_color_to_hex(color_value: str | None) -> str | None:
    if not color_value:
        return None
    parts = [part.strip() for part in color_value.split(",")]
    if len(parts) != 3:
        return None
    try:
        red, green, blue = (int(part) for part in parts)
    except ValueError:
        return None
    if any(channel < 0 or channel > 255 for channel in (red, green, blue)):
        return None
    return f"#{red:02x}{green:02x}{blue:02x}"