"""Tests for UCSC hub parsing."""

from __future__ import annotations

from pathlib import Path

from plotnado import GenomicFigure
from plotnado.hub import UcscHubSession, parse_ucsc_hub
from plotnado.tracks.enums import TrackType


def _write_hub_fixture(tmp_path: Path) -> Path:
    hub_path = tmp_path / "hub.txt"
    genomes_path = tmp_path / "genomes.txt"
    track_db_dir = tmp_path / "grch38"
    track_db_dir.mkdir()
    track_db_path = track_db_dir / "tracksDb.txt"

    hub_path.write_text(
        "\n".join(
            [
                "hub ExampleHub",
                "shortLabel Example hub",
                "longLabel Example hub for tests",
                "genomesFile genomes.txt",
            ]
        ),
        encoding="utf-8",
    )
    genomes_path.write_text(
        "\n".join(
            [
                "genome hg38",
                "trackDb grch38/tracksDb.txt",
            ]
        ),
        encoding="utf-8",
    )
    track_db_path.write_text(
        "\n\n".join(
            [
                "\n".join(
                    [
                        "track containerTrack",
                        "shortLabel Container",
                        "longLabel Container track",
                        "compositeTrack on",
                        "visibility full",
                    ]
                ),
                "\n".join(
                    [
                        "    track signalTrack",
                        "    parent containerTrack",
                        "    shortLabel Signal",
                        "    longLabel Signal coverage",
                        "    type bigWig 0 1000",
                        "    bigDataUrl data/signal.bw",
                        "    color 255,0,0",
                        "    visibility full",
                    ]
                ),
                "\n".join(
                    [
                        "track hiddenContainer",
                        "shortLabel Hidden container",
                        "longLabel Hidden container track",
                        "compositeTrack on",
                        "visibility hide",
                    ]
                ),
                "\n".join(
                    [
                        "    track hiddenBed",
                        "    parent hiddenContainer",
                        "    shortLabel Hidden bed",
                        "    longLabel Hidden bed track",
                        "    type bigBed 6",
                        "    bigDataUrl data/hidden.bb",
                    ]
                ),
                "\n".join(
                    [
                        "    track unsupportedTrack",
                        "    shortLabel Unsupported",
                        "    longLabel Unsupported track",
                        "    type bam",
                        "    bigDataUrl data/reads.bam",
                    ]
                ),
            ]
        ),
        encoding="utf-8",
    )
    return hub_path


def test_parse_ucsc_hub_returns_session(tmp_path: Path) -> None:
    hub_path = _write_hub_fixture(tmp_path)

    result = parse_ucsc_hub(hub_path)

    assert isinstance(result, UcscHubSession)
    assert result.genome == "hg38"
    assert result.short_label == "Example hub"


def test_parse_ucsc_hub_flattens_supported_tracks(tmp_path: Path) -> None:
    hub_path = _write_hub_fixture(tmp_path)

    result = parse_ucsc_hub(hub_path)

    assert [track.name for track in result.template.tracks] == ["signalTrack", "hiddenBed"]

    signal_track, hidden_track = result.template.tracks
    assert signal_track.type == TrackType.BIGWIG
    assert signal_track.path == str((tmp_path / "grch38" / "data" / "signal.bw").resolve())
    assert signal_track.color == "#ff0000"

    assert hidden_track.type == TrackType.BED
    assert hidden_track.path == str((tmp_path / "grch38" / "data" / "hidden.bb").resolve())
    assert hidden_track.height == 0.0


def test_parse_ucsc_hub_can_drop_hidden_tracks(tmp_path: Path) -> None:
    hub_path = _write_hub_fixture(tmp_path)

    result = parse_ucsc_hub(hub_path, include_hidden=False)

    assert [track.name for track in result.template.tracks] == ["signalTrack"]


def test_from_ucsc_hub_builds_figure(tmp_path: Path) -> None:
    hub_path = _write_hub_fixture(tmp_path)

    fig = GenomicFigure.from_ucsc_hub(hub_path)

    assert len(fig.tracks) == 4
    assert fig.tracks[0].__class__.__name__ == "ScaleBar"
    assert fig.tracks[1].__class__.__name__ == "GenomicAxis"
    assert fig.tracks[2].title == "Signal coverage"
    assert fig.tracks[2].height == 1.0
    assert fig.tracks[3].title == "Hidden bed track"
    assert fig.tracks[3].height == 0.0