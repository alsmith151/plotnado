"""Theme primitives for figure-level default styling."""

from enum import Enum

from pydantic import BaseModel, Field, model_validator

from .tracks.base import LabelConfig


class BuiltinTheme(str, Enum):
    """Builtin theme names supported by Plotnado."""

    DEFAULT = "default"
    MINIMAL = "minimal"
    PUBLICATION = "publication"


class Theme(BaseModel):
    """Figure-level theme defaults for tracks and labels."""

    color: str | None = Field(default=None, description="Default color applied to compatible track aesthetics.")
    alpha: float | None = Field(default=None, description="Default opacity applied to compatible track aesthetics.")
    linewidth: float | None = Field(default=None, description="Default line width applied to compatible track aesthetics.")
    font_size: int | None = Field(default=None, description="Default font size for compatible track text elements.")
    font_family: str | None = Field(default=None, description="Global font family applied to labels and text when unset.")
    cmap: str | None = Field(default=None, description="Default colormap for compatible matrix/continuous tracks.")

    highlight_color: str = Field(default="#ffd700", description="Default color used for region highlights.")
    highlight_alpha: float = Field(default=0.15, description="Default opacity used for region highlights.")

    label: LabelConfig = Field(
        default_factory=LabelConfig,
        description="Default label styling configuration merged into tracks.",
    )

    @model_validator(mode="after")
    def _sync_label_fonts(self) -> "Theme":
        if not self.font_family:
            return self

        default_label = LabelConfig()
        if self.label.title_font == default_label.title_font:
            self.label.title_font = self.font_family
        if self.label.scale_font == default_label.scale_font:
            self.label.scale_font = self.font_family
        return self

    @classmethod
    def default(cls) -> "Theme":
        return cls()

    @classmethod
    def minimal(cls) -> "Theme":
        return cls(
            color="#333333",
            alpha=0.8,
            linewidth=1.0,
            font_size=8,
            font_family="DejaVu Sans",
            highlight_color="#cccccc",
            highlight_alpha=0.12,
        )

    @classmethod
    def publication(cls) -> "Theme":
        return cls(
            color="#1f2937",
            alpha=0.9,
            linewidth=1.2,
            font_size=9,
            font_family="DejaVu Sans",
            cmap="viridis",
            highlight_color="#fde68a",
            highlight_alpha=0.2,
            label=LabelConfig(
                title_weight="bold",
                title_size=10,
                scale_size=9,
                label_box_enabled=True,
                label_box_alpha=0.85,
            ),
        )

    @classmethod
    def from_builtin(cls, value: BuiltinTheme | str) -> "Theme":
        """Create a `Theme` from a builtin theme name."""
        builtin = value if isinstance(value, BuiltinTheme) else BuiltinTheme(str(value).lower())

        if builtin == BuiltinTheme.DEFAULT:
            return cls.default()
        if builtin == BuiltinTheme.MINIMAL:
            return cls.minimal()
        if builtin == BuiltinTheme.PUBLICATION:
            return cls.publication()
        raise ValueError(f"Unknown builtin theme: {value}")
