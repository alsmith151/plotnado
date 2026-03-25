"""
Annotation tracks for horizontal and vertical reference lines.
"""

from typing import Any

import matplotlib.axes
from pydantic import Field

from .region import GenomicRegion
from .base import Track
from .enums import TrackType
from .aesthetics import BaseAesthetics
from .registry import registry


class AnnotationAesthetics(BaseAesthetics):
    """
    Defines visual properties for annotation reference lines.

    Inherits color, alpha, and linewidth from BaseAesthetics.

    Attributes:
        linestyle: Matplotlib line style pattern for the line
        zorder: Matplotlib z-order for rendering order
    """

    linestyle: str = Field(default="--", description="Matplotlib line style pattern for the line.")
    zorder: int = Field(default=10, description="Matplotlib z-order used when drawing the line.")


@registry.register(TrackType.HLINE)
class HLineTrack(Track):
    """
    Track for drawing a horizontal reference line.

    Attributes:
        y_value: Y-coordinate for the line
        aesthetics: Visual configuration
    """

    y_value: float = Field(description="Y-axis value where the horizontal line is drawn.")
    aesthetics: AnnotationAesthetics = Field(
        default_factory=AnnotationAesthetics,
        description="Visual styling options for the horizontal reference line.",
    )
    height: float = Field(default=0.0, description="Zero-height overlay track.")

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot horizontal line."""
        ax.axhline(
            y=self.y_value,
            color=self.color,
            linestyle=self.linestyle,
            linewidth=self.linewidth,
            alpha=self.alpha,
            zorder=self.zorder,
        )


@registry.register(TrackType.VLINE)
class VLineTrack(Track):
    """
    Track for drawing a vertical reference line at a genomic position.

    Attributes:
        x_position: Genomic coordinate for the line
        aesthetics: Visual configuration
    """

    x_position: int | str = Field(description="Genomic position where the vertical line is drawn.")
    aesthetics: AnnotationAesthetics = Field(
        default_factory=AnnotationAesthetics,
        description="Visual styling options for the vertical reference line.",
    )
    height: float = Field(default=0.0, description="Zero-height overlay track.")

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot vertical line."""
        from .utils import parse_genomic_value

        x = parse_genomic_value(self.x_position)
        if gr.start <= x <= gr.end:
            ax.axvline(
                x=x,
                color=self.color,
                linestyle=self.linestyle,
                linewidth=self.linewidth,
                alpha=self.alpha,
                zorder=self.zorder,
            )

    def plot_on_axes(self, gr: GenomicRegion, axes: list) -> None:
        """Plot vertical line on multiple axes (global marker)."""
        for ax in axes:
            self.plot(ax, gr)
