"""Theme primitives for figure-level default styling."""

from pydantic import BaseModel, model_validator

from .tracks.base import LabelConfig


class Theme(BaseModel):
    """Figure-level theme defaults for tracks and labels."""

    color: str | None = None
    alpha: float | None = None
    linewidth: float | None = None
    font_size: int | None = None
    font_family: str | None = None
    cmap: str | None = None

    highlight_color: str = "#ffd700"
    highlight_alpha: float = 0.15

    label: LabelConfig = LabelConfig()

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
