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
    palette: list[str] | None = Field(
        default=None,
        description="Optional default palette used for automatic multi-track signal coloring.",
    )
    auto_palette: bool = Field(
        default=False,
        description="Automatically color eligible signal tracks using the theme palette.",
    )
    separator_color: str | None = Field(
        default=None,
        description="Optional horizontal separator color drawn between stacked tracks.",
    )
    separator_alpha: float = Field(default=0.25, description="Opacity of figure-level track separators.")
    separator_linewidth: float = Field(default=0.6, description="Line width of figure-level track separators.")
    subplot_hspace: float = Field(default=0.05, description="Vertical spacing between stacked track axes.")
    margin_left: float = Field(default=0.08, description="Figure subplot left margin fraction.")
    margin_right: float = Field(default=0.985, description="Figure subplot right margin fraction.")
    margin_top: float = Field(default=0.97, description="Figure subplot top margin fraction.")
    margin_bottom: float = Field(default=0.04, description="Figure subplot bottom margin fraction.")

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
            color=None,
            alpha=None,
            linewidth=None,
            font_size=None,
            font_family="DejaVu Sans",
            cmap="magma",
            palette=[
                "#7A3E8E",
                "#2C7FB8",
                "#E15759",
                "#111111",
                "#2CA02C",
                "#3B78B4",
                "#E7298A",
                "#F28E2B",
                "#56B4E9",
                "#0B775E",
            ],
            auto_palette=True,
            separator_color="#b8bec8",
            separator_alpha=0.35,
            separator_linewidth=0.6,
            subplot_hspace=0.03,
            margin_left=0.07,
            margin_right=0.985,
            margin_top=0.975,
            margin_bottom=0.035,
            highlight_color="#e9b2a8",
            highlight_alpha=0.22,
            label=LabelConfig(
                title_weight="normal",
                title_size=11,
                scale_size=11,
                scale_color="#2f2f2f",
                label_box_enabled=False,
                label_box_alpha=0.0,
                title_height=0.86,
                scale_height=0.86,
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
