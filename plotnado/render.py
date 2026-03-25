"""
Render plan and compiler for converting templates to plottable figures.

This module handles the conversion from user-facing templates to internal
rendering structures that can be used by GenomicFigure.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from plotnado.template import Template, TrackSpec
from plotnado.tracks import GenomicRegion
from plotnado.tracks.enums import TrackType


@dataclass
class ResolvedTrack:
    """Resolved template track data ready for figure construction.

    Attributes:
        track_spec: Original template-level track specification.
        index: Zero-based track index in the source template.
        actual_path: Resolved data path after any template compilation step.
        actual_title: Title shown in the figure.
        actual_group: Autoscale group name, if any.
        color_group: Shared color group name, if any.

    Example:
        >>> resolved = ResolvedTrack(track_spec=spec, index=0, actual_title="RNA")
        >>> resolved.to_figure_kwargs()["title"]
        'RNA'
    """

    track_spec: TrackSpec
    index: int

    # Resolved values
    actual_path: Optional[str] = None  # Resolved file path from data source
    actual_title: Optional[str] = None  # Final title to display
    actual_group: Optional[str] = None  # Final group assignment (for autoscale)
    color_group: Optional[str] = None   # Group for color sharing (None if autocolor=False)

    _SOURCE_KWARG_MAP = {
        TrackType.COOLER: "file",
        TrackType.CAPCRUNCHER: "file",
    }

    def get_data(self) -> Optional[str]:
        """Return the resolved data source for the track.

        Returns:
            The configured path or other external data locator, if present.
        """
        return self.track_spec.path or self.actual_path

    def source_kwarg_name(self) -> str:
        """Return the figure-constructor keyword for this track's external data."""
        return self._SOURCE_KWARG_MAP.get(self.track_spec.type, "data")

    def to_figure_kwargs(self) -> dict[str, Any]:
        """Convert the resolved track into ``GenomicFigure`` keyword arguments.

        Returns:
            A flat kwargs dict containing track, aesthetics, and label fields.
            The actual data payload/path is excluded and should be retrieved with
            :meth:`get_data`.

        Example:
            >>> resolved.to_figure_kwargs()["title"]
            'RNA'
        """
        kwargs = {}

        if self.actual_title:
            kwargs['title'] = self.actual_title

        if self.track_spec.style:
            kwargs['style'] = self.track_spec.style

        if self.track_spec.color:
            kwargs['color'] = self.track_spec.color

        if self.track_spec.height != 1.0:
            kwargs['height'] = self.track_spec.height

        if self.actual_group:
            kwargs['autoscale_group'] = self.actual_group

        if self.color_group:
            kwargs['color_group'] = self.color_group

        # Add any additional options
        kwargs.update(self.track_spec.options)

        return kwargs


@dataclass
class RenderPlan:
    """Compiled figure-construction plan derived from a template.

    Attributes:
        template: Original template object.
        tracks: Resolved tracks in plotting order.
        resolved_group_indices: Track references expanded to concrete indices.
        width: Figure width in inches.
        track_height: Default per-track height multiplier.
        genome: Optional genome identifier carried through from the template.
        add_genes: Whether to inject a genes guide track.
        add_axis: Whether to inject a genomic axis guide track.
        add_scalebar: Whether to inject a scale bar guide track.
    """

    template: Template
    tracks: list[ResolvedTrack] = field(default_factory=list)

    # Resolved group indices (title/name -> index mapping), stored here to avoid
    # mutating the source Template object.
    resolved_group_indices: dict[str, list[int]] = field(default_factory=dict)

    # Figure-level settings
    width: float = 12.0
    track_height: float = 1.0
    genome: Optional[str] = None

    # Guide settings
    add_genes: bool = False
    add_axis: bool = True
    add_scalebar: bool = True



class TemplateCompiler:
    """Compiles templates into render plans."""

    @staticmethod
    def compile(template: Template, region: Optional[GenomicRegion] = None) -> RenderPlan:
        """Compile a template into a region-independent render plan.

        Does not mutate ``template``. Resolved group indices are stored in
        :attr:`RenderPlan.resolved_group_indices`.

        Args:
            template: Template specification to compile.
            region: Optional genomic region reserved for future validation hooks.

        Returns:
            A ``RenderPlan`` ready to instantiate a ``GenomicFigure``.

        Example:
            >>> plan = TemplateCompiler.compile(template)
            >>> len(plan.tracks) == len(template.tracks)
            True
        """
        plan = RenderPlan(
            template=template,
            width=template.width,
            track_height=template.track_height,
            genome=template.genome,
            add_genes=template.guides.genes,
            add_axis=template.guides.axis,
            add_scalebar=template.guides.scalebar,
        )

        # Build mapping from track titles/names to indices (case-insensitive)
        title_to_index: dict[str, int] = {}
        for i, track_spec in enumerate(template.tracks):
            key = track_spec.name or track_spec.title or f"Track {i + 1}"
            title_to_index[key.lower()] = i

        # Build group name -> GroupSpec lookup for flag overrides
        group_spec_map = {g.name: g for g in template.groups}

        # Resolve and add tracks
        for i, track_spec in enumerate(template.tracks):
            # Determine autocolor: look up GroupSpec if the track has a group
            autocolor = True
            if track_spec.group and track_spec.group in group_spec_map:
                autocolor = group_spec_map[track_spec.group].autocolor

            resolved = TemplateCompiler._resolve_track(
                track_spec, i, autocolor=autocolor
            )
            plan.tracks.append(resolved)

        # Resolve group track references to indices (stored in plan, not mutating template)
        for group in template.groups:
            resolved_indices = []
            for track_ref in group.tracks:
                ref_lower = track_ref.lower()
                if ref_lower in title_to_index:
                    resolved_indices.append(title_to_index[ref_lower])
                elif track_ref.isdigit():
                    resolved_indices.append(int(track_ref))
                else:
                    raise ValueError(
                        f"Group '{group.name}': track reference '{track_ref}' "
                        f"not found. Available tracks: {list(title_to_index.keys())}"
                    )
            plan.resolved_group_indices[group.name] = resolved_indices

        return plan

    @staticmethod
    def _resolve_track(
        spec: TrackSpec, index: int, autocolor: bool = True
    ) -> "ResolvedTrack":
        """Resolve a track specification to a ready-to-plot track."""
        resolved = ResolvedTrack(
            track_spec=spec,
            index=index,
            actual_path=spec.path,
            actual_title=spec.title or f"Track {index + 1}",
            actual_group=spec.group,
            color_group=spec.group if autocolor else None,
        )

        return resolved
