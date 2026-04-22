"""
Base track classes.
"""

from abc import ABC
from enum import Enum
from typing import Any, get_args, get_origin

import matplotlib.axes
from pydantic import BaseModel, ConfigDict, Field
from pydantic_core import PydanticUndefined

from .enums import DataRangeStyle, FontWeight, Position
from .region import GenomicRegion


class LabelConfig(BaseModel):
    """Unified label configuration for all tracks."""

    plot_title: bool = Field(default=True, description="Show the track title text.")
    plot_scale: bool = Field(default=True, description="Show the y-scale range annotation.")
    label_on_track: bool = Field(
        default=False,
        description="Draw labels inside the track plotting area instead of margins.",
    )
    data_range_style: DataRangeStyle = Field(
        default=DataRangeStyle.TEXT,
        description="How to render the data range annotation for the track.",
    )
    label_box_enabled: bool = Field(
        default=True,
        description="Draw a translucent white box behind label text for legibility.",
    )
    label_box_alpha: float = Field(
        default=0.9,
        description="Opacity of the label background box (0-1).",
    )

    title_location: Position = Field(
        default=Position.LEFT,
        description="Side of the track where the title is anchored.",
    )
    title_height: float = Field(
        default=0.8,
        description="Relative vertical position for title text in data coordinates.",
    )
    title_size: int = Field(default=10, description="Font size for the title label.")
    title_color: str = Field(default="#333333", description="Text color for title label.")
    title_font: str = Field(default="DejaVu Sans", description="Font family for title label.")
    title_weight: FontWeight = Field(
        default=FontWeight.BOLD,
        description="Font weight for title label.",
    )

    scale_location: Position = Field(
        default=Position.RIGHT,
        description="Side of the track where scale text is anchored.",
    )
    scale_height: float = Field(
        default=0.8,
        description="Relative vertical position for scale text in data coordinates.",
    )
    scale_precision: int = Field(
        default=2,
        description="Decimal precision used when formatting scale values.",
    )
    scale_size: int = Field(default=9, description="Font size for scale annotation.")
    scale_color: str = Field(default="#666666", description="Text color for scale annotation.")
    scale_font: str = Field(default="DejaVu Sans", description="Font family for scale annotation.")
    scale_weight: FontWeight = Field(
        default=FontWeight.NORMAL,
        description="Font weight for scale annotation.",
    )

    model_config = ConfigDict(use_enum_values=True)


class Track(BaseModel, ABC):
    """
    Abstract base class for all track types.

    Attributes:
        title: Optional title for the track
        height: Relative height of the track (default 1.0)
    """

    title: str | None = Field(
        default=None,
        description="Display title for the track panel.",
    )
    data: Any | None = Field(
        default=None,
        description="Primary data source consumed by the track implementation.",
    )
    height: float = Field(
        default=1.0,
        description="Relative vertical height allocated to this track.",
    )
    autoscale_group: str | None = Field(
        default=None,
        description="Tracks sharing the same group id are y-scaled together.",
    )
    color_group: str | None = Field(
        default=None,
        description="Tracks sharing this id use the same autocolor/theme-palette color.",
    )
    label: LabelConfig = Field(
        default_factory=LabelConfig,
        description="Unified label and scale rendering configuration.",
    )

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    @classmethod
    def aesthetics_model(cls) -> type[BaseModel] | None:
        field_info = cls.model_fields.get("aesthetics")
        annotation = getattr(field_info, "annotation", None)
        if isinstance(annotation, type) and issubclass(annotation, BaseModel):
            return annotation
        return None

    @classmethod
    def aesthetics_field_names(cls) -> set[str]:
        model = cls.aesthetics_model()
        if model is None:
            return set()
        return set(model.model_fields.keys())

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            return BaseModel.__getattr__(self, name)

        aesthetics = self.__dict__.get("aesthetics")
        if isinstance(aesthetics, BaseModel) and name in aesthetics.__class__.model_fields:
            return getattr(aesthetics, name)
        raise AttributeError(f"{self.__class__.__name__!s} has no attribute {name!r}")

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_"):
            BaseModel.__setattr__(self, name, value)
            return

        if name in type(self).model_fields:
            super().__setattr__(name, value)
            return

        aesthetics = self.__dict__.get("aesthetics")
        if isinstance(aesthetics, BaseModel) and name in aesthetics.__class__.model_fields:
            setattr(aesthetics, name, value)
            return

        super().__setattr__(name, value)

    def has_aesthetic(self, name: str) -> bool:
        aesthetics = getattr(self, "aesthetics", None)
        return isinstance(aesthetics, BaseModel) and name in aesthetics.__class__.model_fields

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
        if isinstance(value, Enum):
            return value.value
        if isinstance(value, (str, int, float, bool, list, dict)):
            return value
        if isinstance(value, BaseModel):
            return value.__class__.__name__
        return value.__class__.__name__

    @staticmethod
    def _enum_choices(annotation: Any) -> list[Any] | None:
        if isinstance(annotation, type) and issubclass(annotation, Enum):
            return [member.value for member in annotation]

        origin = get_origin(annotation)
        if origin is None:
            return None

        choices: list[Any] = []
        for arg in get_args(annotation):
            arg_choices = Track._enum_choices(arg)
            if not arg_choices:
                continue
            for choice in arg_choices:
                if choice not in choices:
                    choices.append(choice)

        return choices or None

    @staticmethod
    def _render_choices(choices: list[Any] | None) -> str:
        if not choices:
            return "—"
        return ", ".join(str(choice) for choice in choices)

    @staticmethod
    def _markdown_cell(value: Any) -> str:
        return str(value).replace("|", "\\|")

    @classmethod
    def options(cls) -> dict[str, dict[str, dict[str, Any]]]:
        """Return programmatically generated option metadata for this track."""
        track_options: dict[str, dict[str, Any]] = {}
        aesthetics_options: dict[str, dict[str, Any]] = {}
        label_options: dict[str, dict[str, Any]] = {}

        for field_name, field_info in cls.model_fields.items():
            if field_name in {"aesthetics", "label"}:
                continue
            track_options[field_name] = {
                "type": cls._render_type_name(field_info.annotation),
                "default": cls._render_default(field_info.default),
                "required": field_info.is_required(),
                "description": field_info.description,
                "choices": cls._enum_choices(field_info.annotation),
            }

        aesthetics_model = cls.aesthetics_model()
        if aesthetics_model is not None:
            for field_name, field_info in aesthetics_model.model_fields.items():
                aesthetics_options[field_name] = {
                    "type": cls._render_type_name(field_info.annotation),
                    "default": cls._render_default(field_info.default),
                    "required": field_info.is_required(),
                    "description": field_info.description,
                    "choices": cls._enum_choices(field_info.annotation),
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
                    "choices": cls._enum_choices(field_info.annotation),
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
            lines.append("| Name | Type | Default | Choices | Required | Description |")
            lines.append("|---|---|---|---|---|---|")
            for name, meta in options[section].items():
                type_cell = cls._markdown_cell(meta["type"])
                default_cell = cls._markdown_cell(meta["default"])
                choices_cell = cls._markdown_cell(cls._render_choices(meta.get("choices")))
                required_cell = cls._markdown_cell(meta["required"])
                description_cell = cls._markdown_cell(meta.get("description") or "—")
                lines.append(
                    f"| {name} | {type_cell} | {default_cell} | {choices_cell} | {required_cell} | {description_cell} |"
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
    data_range_style: DataRangeStyle = DataRangeStyle.TEXT
    label_box_enabled: bool = True
    label_box_alpha: float = 0.9

    title: str = ""
    title_size: int = 10
    title_color: str = "#333333"  # Darker gray
    title_font: str = "DejaVu Sans"
    title_weight: FontWeight = FontWeight.BOLD
    title_location: Position = Position.LEFT
    title_height: float = 0.8

    scale_min: float = 0
    scale_max: float = 1
    scale_precision: int = 2
    scale_size: int = 9
    scale_color: str = "#666666"  # Gray
    scale_font: str = "DejaVu Sans"
    scale_weight: FontWeight = FontWeight.NORMAL
    scale_location: Position = Position.RIGHT
    scale_height: float = 0.8

    model_config = ConfigDict(arbitrary_types_allowed=True, use_enum_values=True)

    @classmethod
    def from_config(
        cls,
        label: LabelConfig,
        gr: GenomicRegion,
        y_min: float,
        y_max: float,
        title: str = "",
        title_color: str | None = None,
    ) -> "TrackLabeller":
        resolved_title_color = label.title_color
        default_label = LabelConfig()
        label_fields_set = getattr(label, "model_fields_set", set())
        if (
            title_color is not None
            and "title_color" not in label_fields_set
            and label.title_color == default_label.title_color
        ):
            resolved_title_color = title_color

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
            title_color=resolved_title_color,
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
            y_pos = self.y_min + (self.y_delta * self.title_height)
            h_align = "left" if self.title_location == "left" else "right"
        title_bbox = self._text_bbox()

        ax.text(
            x_pos,
            y_pos,
            self.title,
            horizontalalignment=h_align,
            verticalalignment="top",
            bbox=title_bbox,
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
        location: Position,
    ) -> None:
        y_min = self._format_scale(self.y_min)
        y_max = self._format_scale(self.y_max)

        x_pos = (
            gr.end - (0.01 * gr.length)
            if location == "right"
            else gr.start + (0.01 * gr.length)
        )
        h_align = "right" if location == "right" else "left"
        scale_bbox = self._text_bbox()

        ax.text(
            x_pos,
            self.y_min + (self.y_delta * self.scale_height),
            f"[ {y_min} - {y_max} ]",
            horizontalalignment=h_align,
            verticalalignment="top",
            bbox=scale_bbox,
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

        title_side: Position
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
