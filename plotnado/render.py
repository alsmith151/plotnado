"""
Render plan and compiler for converting templates to plottable figures.

This module handles the conversion from user-facing templates to internal
rendering structures that can be used by GenomicFigure.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from plotnado.template import Template, TrackSpec, TemplateTrackType
from plotnado.tracks import GenomicRegion


@dataclass
class ResolvedTrack:
    """
    A track specification resolved and ready for plotting.

    Contains all information needed to create a Track object in GenomicFigure.
    """

    track_spec: TrackSpec
    index: int

    # Resolved values
    actual_path: Optional[str] = None  # Resolved file path from data source
    actual_title: Optional[str] = None  # Final title to display
    actual_group: Optional[str] = None  # Final group assignment (for autoscale)
    color_group: Optional[str] = None   # Group for color sharing (None if autocolor=False)

    def get_data(self) -> Optional[str]:
        """Get the data source (file path, URL, or dataframe)."""
        return self.track_spec.path or self.actual_path

    def to_figure_kwargs(self) -> dict[str, Any]:
        """
        Convert to kwargs suitable for GenomicFigure track methods.

        Returns a dict with keys like 'title', 'style', 'color', etc.
        (data is handled separately via get_data())
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
    """
    A compiled plan for rendering a figure.

    Contains all information needed to construct a GenomicFigure and render it
    for a specific region.
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

    def get_track_by_method(self, index: int) -> tuple[str, Any, dict[str, Any]]:
        """
        Get the figure method name, data, and kwargs for a resolved track.

        Returns:
            Tuple of (method_name, data, kwargs) suitable for GenomicFigure
        """
        if index >= len(self.tracks):
            raise IndexError(f"Track index {index} out of range")

        resolved = self.tracks[index]
        spec = resolved.track_spec

        # Map track type to GenomicFigure method (keys are strings due to use_enum_values)
        method_map = {
            TemplateTrackType.BIGWIG.value: 'bigwig',
            TemplateTrackType.BED.value: 'bed',
            TemplateTrackType.NARROWPEAK.value: 'narrowpeak',
            TemplateTrackType.BEDGRAPH.value: 'bigwig',  # bedgraph renders via BigWigTrack; no separate fig.bedgraph() method exists
            TemplateTrackType.GENE.value: 'genes',
            TemplateTrackType.LINKS.value: 'links',
            TemplateTrackType.ANNOTATION.value: 'bed',  # annotation tracks are BED intervals; rendered via fig.bed()
            TemplateTrackType.OVERLAY.value: 'overlay',
            TemplateTrackType.UNKNOWN.value: 'bed',  # Default fallback
        }

        method = method_map.get(str(spec.type), 'bed')
        data = resolved.get_data()
        kwargs = resolved.to_figure_kwargs()

        return method, data, kwargs


class TemplateCompiler:
    """Compiles templates into render plans."""

    @staticmethod
    def compile(template: Template, region: Optional[GenomicRegion] = None) -> RenderPlan:
        """
        Compile a template into a render plan.

        Does NOT mutate the template. Resolved group indices are stored in
        RenderPlan.resolved_group_indices.

        Args:
            template: The template to compile
            region: Optional region for validation

        Returns:
            A RenderPlan ready for rendering
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
