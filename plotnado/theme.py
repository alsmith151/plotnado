"""Theme primitives for figure-level default styling."""

from pydantic import BaseModel

from .tracks.base import LabelConfig


class Theme(BaseModel):
    """Figure-level theme defaults for tracks and labels."""

    color: str | None = None
    alpha: float | None = None
    linewidth: float | None = None
    font_size: int | None = None
    cmap: str | None = None

    highlight_color: str = "#ffd700"
    highlight_alpha: float = 0.15

    label: LabelConfig = LabelConfig()

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
