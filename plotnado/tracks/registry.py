"""
Track registry — single source of truth for track type → class mapping.

All concrete track classes register themselves here using the
``@registry.register`` decorator. This replaces the ``_alias_map`` dict
previously defined in ``GenomicFigure``.

The module-level ``registry`` singleton is populated when each track module
is imported (triggered by ``plotnado/tracks/__init__.py``).

Example:
    Registering a track::

        from plotnado.tracks.registry import registry
        from plotnado.tracks.enums import TrackType

        @registry.register(TrackType.BIGWIG, aliases=["bw", "signal", "bedgraph"])
        class BigWigTrack(Track):
            ...

    Resolving a class::

        entry = registry.get("bigwig")  # or "bw", "signal", "bedgraph"
        track = entry.cls(data="signal.bw", title="H3K27ac")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from .enums import TrackType

if TYPE_CHECKING:
    from .base import Track


@dataclass
class TrackEntry:
    """Metadata for a registered track class.

    Attributes:
        track_type: Canonical ``TrackType`` value for this class.
        cls: The concrete ``Track`` subclass.
        aliases: Additional string aliases accepted by the registry (e.g.,
            ``"bw"`` and ``"signal"`` both resolve to ``BigWigTrack``).
    """

    track_type: TrackType
    cls: type[Track]
    aliases: list[str] = field(default_factory=list)


class TrackRegistry:
    """Registry mapping track type strings to their implementation classes.

    Use the module-level :data:`registry` singleton rather than creating
    instances directly. Track classes register themselves via the
    :meth:`register` decorator at import time.

    Example:
        >>> from plotnado.tracks.registry import registry
        >>> entry = registry.get("bigwig")
        >>> entry.cls.__name__
        'BigWigTrack'
    """

    def __init__(self) -> None:
        self._tracks: dict[str, TrackEntry] = {}

    def register(self, track_type: TrackType, aliases: list[str] = ()):
        """Class decorator that registers a ``Track`` subclass.

        Args:
            track_type: The canonical ``TrackType`` for this class.
            aliases: Additional string keys that also resolve to this class.
                Useful for bedgraph → BigWigTrack, annotation → BedTrack, etc.

        Returns:
            The unmodified class (decorator pattern).

        Example:
            >>> @registry.register(TrackType.BIGWIG, aliases=["bw", "bedgraph"])
            ... class BigWigTrack(Track): ...
        """
        def decorator(cls: type[Track]) -> type[Track]:
            entry = TrackEntry(
                track_type=track_type,
                cls=cls,
                aliases=list(aliases),
            )
            self._tracks[track_type.value] = entry
            for alias in aliases:
                self._tracks[alias] = entry
            return cls
        return decorator

    def get(self, name: str) -> TrackEntry:
        """Retrieve a ``TrackEntry`` by ``TrackType`` value or alias.

        Args:
            name: A ``TrackType`` value string (e.g. ``"bigwig"``) or a
                registered alias (e.g. ``"bw"``).

        Returns:
            The matching ``TrackEntry``.

        Raises:
            KeyError: If ``name`` is not registered.

        Example:
            >>> entry = registry.get("bigwig")
            >>> entry.cls.__name__
            'BigWigTrack'
        """
        if name not in self._tracks:
            available = sorted(self._tracks)
            raise KeyError(
                f"Unknown track type: {name!r}. "
                f"Available: {available}"
            )
        return self._tracks[name]

    def all_types(self) -> list[TrackType]:
        """Return the canonical ``TrackType`` for each registered class.

        Aliases are excluded — each class appears once.

        Returns:
            List of ``TrackType`` values in registration order.
        """
        seen: set[TrackType] = set()
        result: list[TrackType] = []
        for entry in self._tracks.values():
            if entry.track_type not in seen:
                seen.add(entry.track_type)
                result.append(entry.track_type)
        return result

    def all_entries(self) -> dict[str, TrackEntry]:
        """Return the full internal dict (includes aliases).

        Returns:
            Dict mapping every registered string key to its ``TrackEntry``.
        """
        return dict(self._tracks)


#: Module-level singleton populated by ``@registry.register`` decorators
#: when track modules are imported via ``plotnado/tracks/__init__.py``.
registry = TrackRegistry()
