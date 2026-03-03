"""
Base track classes.
"""

from abc import ABC
from typing import Any, Literal

import matplotlib.axes
from pydantic import BaseModel, ConfigDict, Field, model_validator
from pydantic_core import PydanticUndefined

from ._meta import TrackMeta
from .region import GenomicRegion


class LabelConfig(BaseModel):
    """Unified label configuration for all tracks."""

    plot_title: bool = True
    plot_scale: bool = True
    label_on_track: bool = False
    data_range_style: Literal["text", "colorbar", "none"] = "text"
    label_box_enabled: bool = True
    label_box_alpha: float = 0.9

    title_location: Literal["left", "right"] = "left"
    title_height: float = 0.8
    title_size: int = 10
    title_color: str = "#333333"
    title_font: str = "DejaVu Sans"
    title_weight: Literal["normal", "bold"] = "bold"

    scale_location: Literal["left", "right"] = "right"
    scale_height: float = 0.8
    scale_precision: int = 2
    scale_size: int = 9
    scale_color: str = "#666666"
    scale_font: str = "DejaVu Sans"
    scale_weight: Literal["normal", "bold"] = "normal"


class Track(BaseModel, ABC, metaclass=TrackMeta):
    """
    Abstract base class for all track types.

    Attributes:
        title: Optional title for the track
        height: Relative height of the track (default 1.0)
    """

    title: str | None = None
    data: Any | None = None
    height: float = 1.0
    autoscale_group: str | None = None
    label: LabelConfig = Field(default_factory=LabelConfig)

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def _merge_nested_aesthetics(cls, data: Any) -> Any:
        """Allow nested aesthetics input while exposing flattened kwargs."""
        if not isinstance(data, dict) or "aesthetics" not in data:
            return data

        aesthetics_data = data.get("aesthetics")
        if isinstance(aesthetics_data, BaseModel):
            aesthetics_data = aesthetics_data.model_dump()

        if not isinstance(aesthetics_data, dict):
            return data

        merged = dict(data)
        for key, value in aesthetics_data.items():
            if key in cls.model_fields and key != "aesthetics":
                merged.setdefault(key, value)
        return merged

    @model_validator(mode="after")
    def _sync_aesthetics_model(self) -> "Track":
        """Keep `self.aesthetics` aligned with flattened aesthetic fields."""
        aesthetics_class = getattr(type(self), "__track_aesthetics_class__", None)
        if aesthetics_class is None or "aesthetics" not in type(self).model_fields:
            return self

        aesthetics_payload: dict[str, Any] = {}
        current_aesthetics = getattr(self, "aesthetics", None)
        if isinstance(current_aesthetics, BaseModel):
            aesthetics_payload.update(current_aesthetics.model_dump())
        elif isinstance(current_aesthetics, dict):
            aesthetics_payload.update(current_aesthetics)

        for field_name in getattr(
            type(self), "__track_flattened_aesthetics_fields__", set()
        ):
            if field_name in type(self).model_fields:
                aesthetics_payload[field_name] = getattr(self, field_name)

        object.__setattr__(self, "aesthetics", aesthetics_class(**aesthetics_payload))
        return self

    def fetch_data(self, gr: GenomicRegion) -> Any:
        """Fetch data for the given genomic region."""
        raise NotImplementedError

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the track on the given axes for the given region."""
        raise NotImplementedError

    @staticmethod
    def _render_type_name(annotation: Any) -> str:
        if annotation is None:
            return "None"
        if isinstance(annotation, type):
            return annotation.__name__
        return str(annotation).replace("typing.", "")

    @staticmethod
    def _render_default(value: Any) -> Any:
        if value is None:
            return "—"
        if value is PydanticUndefined:
            return "*(required)*"
        if isinstance(value, (str, int, float, bool, list, dict)):
            return value
        if isinstance(value, BaseModel):
            return value.__class__.__name__
        return value.__class__.__name__

    @classmethod
    def options(cls) -> dict[str, dict[str, dict[str, Any]]]:
        """Return programmatically generated option metadata for this track."""
        native_fields = set(getattr(cls, "__track_native_fields__", set()))
        flattened_fields = set(
            getattr(cls, "__track_flattened_aesthetics_fields__", set())
        )

        track_options: dict[str, dict[str, Any]] = {}
        aesthetics_options: dict[str, dict[str, Any]] = {}
        label_options: dict[str, dict[str, Any]] = {}

        for field_name, field_info in cls.model_fields.items():
            if field_name in {"aesthetics", "label"}:
                continue
            section = (
                aesthetics_options
                if field_name in flattened_fields and field_name not in native_fields
                else track_options
            )
            section[field_name] = {
                "type": cls._render_type_name(field_info.annotation),
                "default": cls._render_default(field_info.default),
                "required": field_info.is_required(),
                "description": field_info.description,
            }

        label_info = cls.model_fields.get("label")
        if label_info and isinstance(label_info.annotation, type) and issubclass(
            label_info.annotation, BaseModel
        ):
            for field_name, field_info in label_info.annotation.model_fields.items():
                label_options[field_name] = {
                    "type": cls._render_type_name(field_info.annotation),
                    "default": cls._render_default(field_info.default),
                    "required": field_info.is_required(),
                    "description": field_info.description,
                }

        return {
            "track": track_options,
            "aesthetics": aesthetics_options,
            "label": label_options,
        }

    @classmethod
    def options_markdown(cls) -> str:
        """Return options metadata as markdown tables for notebook display."""
        options = cls.options()
        lines = [f"## {cls.__name__} options", ""]

        for section in ("track", "aesthetics", "label"):
            if not options[section]:
                continue
            lines.append(f"### {section.title()} fields")
            lines.append("| Name | Type | Default | Required |")
            lines.append("|---|---|---|---|")
            for name, meta in options[section].items():
                lines.append(
                    f"| {name} | {meta['type']} | {meta['default']} | {meta['required']} |"
                )
            lines.append("")

        return "\n".join(lines)

    @staticmethod
    def _truncate_repr(value: str, max_len: int = 60) -> str:
        if len(value) <= max_len:
            return value
        return f"{value[: max_len - 3]}..."

    def __repr__(self) -> str:
        parts: list[str] = []
        if self.title:
            parts.append(f"title={self.title!r}")
        if self.data is not None:
            data_text = self._truncate_repr(str(self.data), max_len=60)
            parts.append(f"data={data_text!r}")
        body = ", ".join(parts)
        return f"{self.__class__.__name__}({body})"

    def _repr_html_(self) -> str:
        import html

        rows: list[tuple[str, str]] = [("Type", self.__class__.__name__)]
        if self.title:
            rows.append(("Title", str(self.title)))
        if self.data is not None:
            rows.append(("Data", self._truncate_repr(str(self.data), max_len=120)))
        rows.append(("Height", str(self.height)))
        if self.autoscale_group:
            rows.append(("Autoscale group", self.autoscale_group))

        row_html = "".join(
            f"<tr><th style='text-align:left;padding-right:12px'>{html.escape(key)}</th><td>{html.escape(value)}</td></tr>"
            for key, value in rows
        )
        return f"<table>{row_html}</table>"


class TrackLabeller(BaseModel):
    """
    Handles track labelling (title and scale display).
    """

    gr: GenomicRegion
    y_min: float
    y_max: float

    plot_title: bool = True
    plot_scale: bool = True
    label_on_track: bool = False
    data_range_style: Literal["text", "colorbar", "none"] = "text"
    label_box_enabled: bool = True
    label_box_alpha: float = 0.9

    title: str = ""
    title_size: int = 10
    title_color: str = "#333333"  # Darker gray
    title_font: str = "DejaVu Sans"
    title_weight: Literal["normal", "bold"] = "bold"
    title_location: Literal["left", "right"] = "left"
    title_height: float = 0.8

    scale_min: float = 0
    scale_max: float = 1
    scale_precision: int = 2
    scale_size: int = 9
    scale_color: str = "#666666"  # Gray
    scale_font: str = "DejaVu Sans"
    scale_weight: Literal["normal", "bold"] = "normal"
    scale_location: Literal["left", "right"] = "right"
    scale_height: float = 0.8

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @classmethod
    def from_config(
        cls,
        label: LabelConfig,
        gr: GenomicRegion,
        y_min: float,
        y_max: float,
        title: str = "",
    ) -> "TrackLabeller":
        return cls(
            gr=gr,
            y_min=y_min,
            y_max=y_max,
            plot_title=label.plot_title,
            plot_scale=label.plot_scale,
            label_on_track=label.label_on_track,
            data_range_style=label.data_range_style,
            label_box_enabled=label.label_box_enabled,
            label_box_alpha=label.label_box_alpha,
            title=title,
            scale_min=y_min,
            scale_max=y_max,
            title_location=label.title_location,
            title_height=label.title_height,
            title_size=label.title_size,
            title_color=label.title_color,
            title_font=label.title_font,
            title_weight=label.title_weight,
            scale_location=label.scale_location,
            scale_height=label.scale_height,
            scale_precision=label.scale_precision,
            scale_size=label.scale_size,
            scale_color=label.scale_color,
            scale_font=label.scale_font,
            scale_weight=label.scale_weight,
        )

    @staticmethod
    def _default_text_bbox(alpha: float) -> dict:
        return {
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": alpha,
            "boxstyle": "round,pad=0.2",
        }

    def _text_bbox(self) -> dict | None:
        if not self.label_box_enabled:
            return None
        return self._default_text_bbox(self.label_box_alpha)

    @property
    def y_delta(self) -> float:
        return self.y_max - self.y_min

    def _plot_title(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        if self.label_on_track:
            x_pos = gr.start + (0.02 * gr.length)
            y_pos = self.y_min + (self.y_delta * self.title_height)
            h_align = "left"
        else:
            x_pos = (
                gr.start + (0.01 * gr.length)
                if self.title_location == "left"
                else gr.end - (0.01 * gr.length)
            )
            y_pos = self.y_delta * self.title_height
            h_align = "left" if self.title_location == "left" else "right"

        ax.text(
            x_pos,
            y_pos,
            self.title,
            horizontalalignment=h_align,
            verticalalignment="top",
            bbox=self._text_bbox(),
            fontdict={
                "size": self.title_size,
                "color": self.title_color,
                "fontname": self.title_font,
                "weight": self.title_weight,
            },
        )

    def _format_scale(self, value: float) -> str:
        """Format scale value for display."""
        if value % 1 == 0:
            return str(int(value))
        return f"{value:.{self.scale_precision}f}"

    def _plot_scale(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        self._plot_scale_at(ax, gr, self.scale_location)

    def _plot_scale_at(
        self,
        ax: matplotlib.axes.Axes,
        gr: GenomicRegion,
        location: Literal["left", "right"],
    ) -> None:
        y_min = self._format_scale(self.y_min)
        y_max = self._format_scale(self.y_max)

        x_pos = (
            gr.end - (0.01 * gr.length)
            if location == "right"
            else gr.start + (0.01 * gr.length)
        )
        h_align = "right" if location == "right" else "left"

        ax.text(
            x_pos,
            self.y_delta * self.scale_height,
            f"[ {y_min} - {y_max} ]",
            horizontalalignment=h_align,
            verticalalignment="top",
            bbox=self._text_bbox(),
            fontdict={
                "size": self.scale_size,
                "color": self.scale_color,
                "fontname": self.scale_font,
                "weight": self.scale_weight,
            },
        )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> "TrackLabeller":
        should_plot_scale = self.plot_scale and self.data_range_style == "text"
        should_plot_title = self.plot_title and bool(self.title)

        title_side: Literal["left", "right"]
        if self.label_on_track:
            title_side = "left"
        else:
            title_side = self.title_location

        if should_plot_title:
            self._plot_title(ax, gr)
        if should_plot_scale:
            original_scale_location = self.scale_location
            if should_plot_title:
                self.scale_location = "right" if title_side == "left" else "left"
            try:
                self._plot_scale(ax, gr)
            finally:
                self.scale_location = original_scale_location

        import plotnado.tracks as pg_tracks

        pg_tracks.clean_axis(ax)
        return self
