"""
Option helpers for aesthetics discoverability.
"""

from typing import Any

from .base import Track


def list_options(track_cls: type[Track]) -> dict[str, dict[str, dict[str, Any]]]:
    """Return generated option metadata for a Track subclass.

    This is notebook-friendly and avoids guessing kwargs.
    """
    return track_cls.options()
