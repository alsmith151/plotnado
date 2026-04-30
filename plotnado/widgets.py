"""Notebook widgets for interactive figure inspection."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from .tracks import HighlightsFromFile, HLineTrack, Track, VLineTrack

if TYPE_CHECKING:
    from .figure import GenomicFigure
    from .tracks import GenomicRegion


def _require_ipywidgets() -> Any:
    try:
        import ipywidgets as widgets
    except ImportError as exc:  # pragma: no cover - import error path depends on environment
        raise ImportError(
            "Track visibility widgets require ipywidgets. Install plotnado[notebook] "
            "or run `uv pip install ipywidgets`."
        ) from exc
    return widgets


def _is_toggleable_track(track: Track) -> bool:
    return not isinstance(track, (HighlightsFromFile, HLineTrack, VLineTrack))


def _default_restore_height(track: Track) -> float:
    field_info = type(track).model_fields.get("height")
    default_height = getattr(field_info, "default", None)
    if isinstance(default_height, (int, float)) and default_height > 0:
        return float(default_height)
    return 1.0


def _track_label(track: Track) -> str:
    return track.title or track.__class__.__name__


@dataclass
class TrackToggleControl:
    """Pair a track with its notebook checkbox control."""

    track: Track
    label: str
    checkbox: Any


class TrackVisibilityWidget:
    """Interactive notebook widget for showing and hiding figure tracks."""

    def __init__(self, figure: GenomicFigure, region: str | GenomicRegion):
        widgets = _require_ipywidgets()

        self.figure = figure
        self.region = region
        self.controls: list[TrackToggleControl] = []
        self._widgets = widgets
        self._cached_heights: dict[int, float] = {}
        self.output = widgets.Output()

        checkbox_widgets: list[Any] = []
        for track in self.figure.tracks:
            if not _is_toggleable_track(track):
                continue
            self._cached_heights[id(track)] = (
                track.height if track.height > 0 else _default_restore_height(track)
            )
            checkbox = widgets.Checkbox(
                value=track.height > 0,
                description=_track_label(track),
                indent=False,
            )
            checkbox.observe(self._build_toggle_handler(track), names="value")
            self.controls.append(
                TrackToggleControl(track=track, label=checkbox.description, checkbox=checkbox)
            )
            checkbox_widgets.append(checkbox)

        self.widget = widgets.VBox(
            [
                widgets.HTML("<strong>Track visibility</strong>"),
                widgets.VBox(checkbox_widgets),
                self.output,
            ]
        )
        self.refresh()

    def _build_toggle_handler(self, track: Track):
        def _handle(change: dict[str, Any]) -> None:
            visible = bool(change["new"])
            if visible:
                track.height = self._cached_heights.get(id(track), _default_restore_height(track))
            else:
                if track.height > 0:
                    self._cached_heights[id(track)] = track.height
                track.height = 0.0
            self.refresh()

        return _handle

    def refresh(self) -> None:
        from IPython.display import clear_output, display

        with self.output:
            clear_output(wait=True)
            if not any(track.height > 0 for track in self.figure.tracks if _is_toggleable_track(track)):
                display(self._widgets.HTML("<em>No visible panel tracks.</em>"))
                return

            rendered = self.figure.plot(self.region)
            if rendered is None:
                display(self._widgets.HTML("<em>No visible panel tracks.</em>"))
                return
            display(rendered)

    def _ipython_display_(self) -> None:  # pragma: no cover - exercised by notebooks
        from IPython.display import display

        display(self.widget)