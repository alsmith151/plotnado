"""
Shared base aesthetics models and option helpers for track discoverability.
"""

from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from .base import Track


class BaseAesthetics(BaseModel):
    """Shared visual properties for single-color track types.

    All track aesthetics classes that render a single primary color should
    inherit from this base. Provides ``color``, ``alpha``, and ``linewidth``
    as common fields so theme application and introspection have a consistent
    interface.

    Example:
        >>> class MyAesthetics(BaseAesthetics):
        ...     style: str = "fill"
    """

    model_config = ConfigDict(use_enum_values=True)

    color: str = Field(default="steelblue", description="Primary render color.")
    alpha: float = Field(default=1.0, ge=0.0, le=1.0, description="Opacity (0–1).")
    linewidth: float = Field(default=1.0, ge=0.0, description="Stroke width in points.")


class BaseMultiColorAesthetics(BaseModel):
    """Shared visual properties for multi-color track types (collections, overlays).

    Used instead of ``BaseAesthetics`` when a track renders multiple series,
    each needing its own color (e.g., BigWigCollection, OverlayTrack).

    Example:
        >>> class MyOverlayAesthetics(BaseMultiColorAesthetics):
        ...     show_labels: bool = True
    """

    model_config = ConfigDict(use_enum_values=True)

    colors: list[str] | None = Field(
        default=None,
        description="Per-series color list. Falls back to theme palette when None.",
    )
    alpha: float = Field(default=1.0, ge=0.0, le=1.0, description="Opacity (0–1).")
    linewidth: float = Field(default=1.0, ge=0.0, description="Stroke width in points.")


def list_options(track_cls: type[Track]) -> dict[str, dict[str, dict[str, Any]]]:
    """Return generated option metadata for a Track subclass.

    Useful for interactive exploration and avoiding guessed kwargs.

    Args:
        track_cls: A concrete subclass of :class:`~plotnado.tracks.base.Track`.

    Returns:
        Nested dict mapping section (``"track"``, ``"aesthetics"``, ``"label"``)
        to field name to metadata dict.

    Example:
        >>> from plotnado.tracks.bigwig import BigWigTrack
        >>> opts = list_options(BigWigTrack)
        >>> "color" in opts["aesthetics"]
        True
    """
    return track_cls.options()
