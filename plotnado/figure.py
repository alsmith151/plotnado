"""GenomicFigure class for composing and plotting genomic tracks."""

from __future__ import annotations

import importlib.resources
import json
import math
from pathlib import Path
from typing import Any, Self, Unpack, overload

import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger
from pydantic import BaseModel

from .tracks import (
    BedTrack,
    BigWigCollection,
    BigWigDiff,
    BigWigTrack,
    CapcruncherTrack,
    CoolerAverage,
    CoolerTrack,
    Genes,
    GenomicAxis,
    GenomicRegion,
    HighlightsFromFile,
    HLineTrack,
    LinksTrack,
    NarrowPeakTrack,
    OverlayTrack,
    ScaleBar,
    Spacer,
    Track,
    VLineTrack,
    BigwigOverlay,
    LabelConfig,
    list_options,
)
from .tracks.enums import (
    BigWigDiffMethod,
    CollectionStyle,
    CoolerTransform,
    DataRangeStyle,
    DisplayMode,
    FontWeight,
    GeneLabelOverlapStrategy,
    GeneLabelStyle,
    NarrowPeakColorBy,
    PlotStyle,
    Position,
)
from .theme import BuiltinTheme, Theme
from ._kwargs import (
    AxisKwargs,
    BedKwargs,
    BigwigCollectionKwargs,
    BigwigDiffKwargs,
    BigwigKwargs,
    BigwigOverlayKwargs,
    CoolerKwargs,
    GenesKwargs,
    HighlightsKwargs,
    HlineKwargs,
    LinksKwargs,
    NarrowpeakKwargs,
    OverlayKwargs,
    ScalebarKwargs,
    SpacerKwargs,
    VlineKwargs,
)


class GenomicFigure:
    """Compose and plot multiple genomic tracks.

    Example:
        fig = (
            GenomicFigure()
            .bigwig("signal.bw", title="H3K27ac")
            .genes()
            .axis()
            .scalebar()
        )
        fig.plot("chr1:1000000-2000000")
    """

    def __init__(
        self,
        tracks: list[Track] | None = None,
        width: float = 12,
        track_height: float = 2.0,
        theme: Theme | BuiltinTheme | str | None = None,
    ):
        """Create a figure container for genomic tracks.

        Args:
            tracks: Optional pre-populated list of tracks.
            width: Figure width in inches.
            track_height: Height scaling factor per track.
            theme: Optional `Theme` (or builtin name like "default", "minimal",
                "publication") applied to tracks on add/initialization.
        """
        self.tracks: list[Track] = tracks or []
        self.width = width
        self.track_height = track_height
        self.theme = self._resolve_theme(theme)
        self.highlight_color = (
            self.theme.highlight_color if self.theme is not None else "#ffd700"
        )
        self.highlight_alpha = self.theme.highlight_alpha if self.theme is not None else 0.15
        self._highlight_regions: list[GenomicRegion] = []
        self._autocolor_palette: str | None = None
        self._autoscale: bool = False

        if self.theme is not None:
            self.tracks = [self._apply_theme_to_track(track) for track in self.tracks]

    @staticmethod
    def _resolve_theme(theme: Theme | BuiltinTheme | str | None) -> Theme | None:
        if theme is None:
            return None
        if isinstance(theme, Theme):
            return theme
        try:
            return Theme.from_builtin(theme)
        except ValueError as exc:
            raise ValueError(
                f"Unknown builtin theme: {theme}. Available: {[item.value for item in BuiltinTheme]}"
            ) from exc

    def add_track(self, track: str | Track, **kwargs: Any) -> Self:
        """Add a track instance or track alias to the figure.

        Args:
            track: Track instance or alias (for example `"bigwig"`, `"genes"`).
            **kwargs: Parameters used when `track` is provided as an alias.

        Returns:
            Self, enabling method chaining.
        """
        if isinstance(track, str):
            track = self._create_track_from_alias(track, **kwargs)
        track = self._apply_theme_to_track(track)
        track = self._apply_autocolor_to_track(track)
        self.tracks.append(track)
        return self

    def _apply_autocolor_to_track(self, track: Track) -> Track:
        if self._autocolor_palette is None or not self._should_autocolor_track(track):
            return track

        import matplotlib.colors as mcolors

        cmap = plt.get_cmap(self._autocolor_palette)
        index = sum(1 for existing_track in self.tracks if self._should_autocolor_track(existing_track))
        track.color = mcolors.to_hex(cmap(index % cmap.N))
        return track

    @staticmethod
    def _is_meta_track(track: Track) -> bool:
        return isinstance(
            track,
            (
                ScaleBar,
                GenomicAxis,
                HighlightsFromFile,
                Spacer,
                HLineTrack,
                VLineTrack,
            ),
        )

    @classmethod
    def _should_autocolor_track(cls, track: Track) -> bool:
        return track.has_aesthetic("color") and not cls._is_meta_track(track)

    # BEGIN AUTO-GENERATED OVERLOAD: bigwig
    @overload
    def bigwig(
        self,
        data: Any,
        /,
        *,
        autoscale_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        y_max: float | None = ...,
        y_min: float | None = ...,
        alpha: float = ...,
        color: str = ...,
        fill: bool = ...,
        linewidth: float = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        scatter_point_size: float = ...,
        style: PlotStyle = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: bigwig
    def bigwig(self, data: Any, /, **kwargs: Unpack[BigwigKwargs]) -> Self:
        """Add a BigWig signal track.

        Args:
            data: BigWig data source, typically a file path or compatible dataframe.
            **kwargs: `BigWigTrack` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("bigwig", data=data, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: genes
    @overload
    def genes(
        self,
        genome: str = "hg38",
        /,
        *,
        autoscale_group: str | None = ...,
        data: Path | str | pd.DataFrame | None = ...,
        gene_count: int = ...,
        height: float = ...,
        row_scale: float = ...,
        small_relative: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        arrow_size: float = ...,
        chevron_height_ratio: float = ...,
        chevron_margin_bp: float = ...,
        chevron_max_count: int = ...,
        chevron_min_width_bp: float = ...,
        chevron_target_spacing_bp: float = ...,
        chevron_vertical_offset_ratio: float = ...,
        chevron_width_fraction: float = ...,
        color: str = ...,
        display: DisplayMode = ...,
        exon_color: str = ...,
        exon_edge_color: str = ...,
        exon_linewidth: float = ...,
        fill: bool = ...,
        gene_label_font_size: int = ...,
        gene_label_style: GeneLabelStyle = ...,
        interval_height: float = ...,
        intron_color: str = ...,
        intron_linewidth: float = ...,
        label_connector_linewidth: float = ...,
        label_connectors: bool = ...,
        label_max_chars: int = ...,
        label_min_overlap_bp: int = ...,
        label_min_overlap_fraction: float = ...,
        label_offset_fraction: float = ...,
        label_overlap_strategy: GeneLabelOverlapStrategy = ...,
        label_stagger_offset: float = ...,
        label_vertical_offset: float = ...,
        max_number_of_rows: int = ...,
        minimum_gene_length: int = ...,
        show_labels: bool = ...,
        style: PlotStyle = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: genes
    def genes(self, genome: str = "hg38", /, **kwargs: Unpack[GenesKwargs]) -> Self:
        """Add a genes annotation track.

        Args:
            genome: Genome identifier used for bundled gene annotations.
            **kwargs: `Genes` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("genes", genome=genome, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: axis
    @overload
    def axis(
        self,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        show_chromosome: bool = ...,
        title: str = ...,
        axis_linewidth: float = ...,
        chromosome_fontweight: FontWeight = ...,
        color: str = ...,
        font_size: int = ...,
        num_ticks: int = ...,
        tick_color: str = ...,
        tick_height: float = ...,
        tick_linewidth: float = ...,
        use_human_readable_labels: bool = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: axis
    def axis(self, **kwargs: Unpack[AxisKwargs]) -> Self:
        """Add a genomic coordinate axis track.

        Args:
            **kwargs: `GenomicAxis` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("axis", **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: scalebar
    @overload
    def scalebar(
        self,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str = ...,
        bar_linewidth: float = ...,
        color: str = ...,
        font_size: int = ...,
        label_offset: float = ...,
        position: Position = ...,
        scale_distance: float | None = ...,
        style: PlotStyle = ...,
        tick_height: float = ...,
        tick_linewidth: float = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: scalebar
    def scalebar(self, **kwargs: Unpack[ScalebarKwargs]) -> Self:
        """Add a genomic scale bar track.

        Args:
            **kwargs: `ScaleBar` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("scalebar", **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: spacer
    @overload
    def spacer(
        self,
        height: float = 0.5,
        /,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        title: str = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: spacer
    def spacer(self, height: float = 0.5, **kwargs: Unpack[SpacerKwargs]) -> Self:
        """Add an empty spacer track.

        Args:
            height: Spacer height used to separate neighboring tracks.
            **kwargs: `Spacer` constructor kwargs; shorthand composition routes
                label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("spacer", height=height, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: bed
    @overload
    def bed(
        self,
        data: Any,
        /,
        *,
        autoscale_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        color: str = ...,
        display: DisplayMode = ...,
        edge_color: str = ...,
        font_size: int = ...,
        interval_height: float = ...,
        label_field: str = ...,
        max_rows: int = ...,
        rect_linewidth: float = ...,
        show_labels: bool = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: bed
    def bed(self, data: Any, /, **kwargs: Unpack[BedKwargs]) -> Self:
        """Add a BED interval track.

        Args:
            data: BED data source, typically a file path or compatible dataframe.
            **kwargs: `BedTrack` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("bed", data=data, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: cooler
    @overload
    def cooler(
        self,
        file: str,
        /,
        *,
        autoscale_group: str | None = ...,
        balance: bool = ...,
        data: Any | None = ...,
        height: float = ...,
        resolution: int | None = ...,
        title: str | None = ...,
        transform: CoolerTransform = ...,
        cmap: str = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: cooler
    def cooler(self, file: str, /, **kwargs: Unpack[CoolerKwargs]) -> Self:
        """Add a cooler/mcool matrix track.

        Args:
            file: Cooler file path.
            **kwargs: `CoolerTrack` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("cooler", file=file, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: bigwig_collection
    @overload
    def bigwig_collection(
        self,
        files: list[str],
        /,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        colors: list[str] | None = ...,
        labels: list[str] | None = ...,
        style: CollectionStyle = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: bigwig_collection
    def bigwig_collection(
        self, files: list[str], /, **kwargs: Unpack[BigwigCollectionKwargs]
    ) -> Self:
        """Add a collection track for multiple BigWig files.

        Args:
            files: List of BigWig file paths.
            **kwargs: `BigWigCollection` constructor kwargs; shorthand composition
                routes aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("bigwig_collection", files=files, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: bigwig_diff
    @overload
    def bigwig_diff(
        self,
        file_a: str,
        file_b: str,
        /,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        method: BigWigDiffMethod = ...,
        title: str | None = ...,
        bar_alpha: float = ...,
        linewidth: float = ...,
        negative_color: str = ...,
        positive_color: str = ...,
        zero_line_alpha: float = ...,
        zero_line_color: str = ...,
        zero_line_width: float = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: bigwig_diff
    def bigwig_diff(
        self, file_a: str, file_b: str, /, **kwargs: Unpack[BigwigDiffKwargs]
    ) -> Self:
        """Add a two-signal BigWig difference track.

        Args:
            file_a: First BigWig file path.
            file_b: Second BigWig file path.
            **kwargs: `BigWigDiff` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("bigwig_diff", file_a=file_a, file_b=file_b, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: bigwig_overlay
    @overload
    def bigwig_overlay(
        self,
        tracks: list[Any],
        /,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        colors: list[str] | None = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        show_labels: bool = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: bigwig_overlay
    def bigwig_overlay(
        self, tracks: list[Any], /, **kwargs: Unpack[BigwigOverlayKwargs]
    ) -> Self:
        """Add an overlay track containing multiple BigWig signals.

        Args:
            tracks: Track inputs for overlay, usually paths or `BigWigTrack` objects.
            **kwargs: `BigwigOverlay` constructor kwargs; shorthand composition
                routes aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("bigwig_overlay", tracks=tracks, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: overlay
    @overload
    def overlay(
        self,
        tracks: list[Any],
        /,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        colors: list[str] | None = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        show_labels: bool = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: overlay
    def overlay(self, tracks: list[Any], /, **kwargs: Unpack[OverlayKwargs]) -> Self:
        """Add a generic overlay track for multiple tracks on one axis.

        Args:
            tracks: Track inputs for overlay, usually track instances or BigWig paths.
            **kwargs: `OverlayTrack` constructor kwargs; shorthand composition
                routes aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("overlay", tracks=tracks, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: narrowpeak
    @overload
    def narrowpeak(
        self,
        data: Any,
        /,
        *,
        autoscale_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        cmap: str = ...,
        color: str = ...,
        color_by: NarrowPeakColorBy | None = ...,
        display: DisplayMode = ...,
        edge_color: str = ...,
        font_size: int = ...,
        interval_height: float = ...,
        label_field: str = ...,
        max_rows: int = ...,
        max_score: float | None = ...,
        min_score: float | None = ...,
        rect_linewidth: float = ...,
        show_labels: bool = ...,
        show_summit: bool = ...,
        summit_color: str = ...,
        summit_width: float = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: narrowpeak
    def narrowpeak(self, data: Any, /, **kwargs: Unpack[NarrowpeakKwargs]) -> Self:
        """Add a narrowPeak interval track.

        Args:
            data: NarrowPeak data source, typically a file path or compatible dataframe.
            **kwargs: `NarrowPeakTrack` constructor kwargs; shorthand composition
                routes aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("narrowpeak", data=data, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: links
    @overload
    def links(
        self,
        data: Any,
        /,
        *,
        autoscale_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        cmap: str = ...,
        color: str = ...,
        color_by_score: bool = ...,
        edge_color: str | None = ...,
        linewidth: float = ...,
        max_height: float = ...,
        max_score: float | None = ...,
        min_score: float | None = ...,
        y_baseline: float = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: links
    def links(self, data: Any, /, **kwargs: Unpack[LinksKwargs]) -> Self:
        """Add a genomic links/arcs track.

        Args:
            data: BEDPE-like data source for links.
            **kwargs: `LinksTrack` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("links", data=data, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: highlights
    @overload
    def highlights(
        self,
        data: Any,
        /,
        *,
        autoscale_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        color: str = ...,
        edge_color: str | None = ...,
        linewidth: float = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: highlights
    def highlights(self, data: Any, /, **kwargs: Unpack[HighlightsKwargs]) -> Self:
        """Add file-based highlighted regions spanning tracks.

        Args:
            data: Highlight interval data source.
            **kwargs: `HighlightsFromFile` constructor kwargs; shorthand composition
                routes aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("highlight", data=data, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: hline
    @overload
    def hline(
        self,
        y_value: float,
        /,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        color: str = ...,
        linestyle: str = ...,
        linewidth: float = ...,
        zorder: int = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: hline
    def hline(self, y_value: float, /, **kwargs: Unpack[HlineKwargs]) -> Self:
        """Add a horizontal reference line.

        Args:
            y_value: Y value where the line is drawn.
            **kwargs: `HLineTrack` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("hline", y_value=y_value, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: vline
    @overload
    def vline(
        self,
        x_position: int | str,
        /,
        *,
        autoscale_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        color: str = ...,
        linestyle: str = ...,
        linewidth: float = ...,
        zorder: int = ...,
        data_range_style: DataRangeStyle = ...,
        label_box_alpha: float = ...,
        label_box_enabled: bool = ...,
        label_on_track: bool = ...,
        plot_scale: bool = ...,
        plot_title: bool = ...,
        scale_color: str = ...,
        scale_font: str = ...,
        scale_height: float = ...,
        scale_location: Position = ...,
        scale_precision: int = ...,
        scale_size: int = ...,
        scale_weight: FontWeight = ...,
        title_color: str = ...,
        title_font: str = ...,
        title_height: float = ...,
        title_location: Position = ...,
        title_size: int = ...,
        title_weight: FontWeight = ...,
    ) -> Self: ...
    # END AUTO-GENERATED OVERLOAD: vline
    def vline(self, x_position: int | str, /, **kwargs: Unpack[VlineKwargs]) -> Self:
        """Add a vertical reference line.

        Args:
            x_position: Genomic coordinate for the line.
            **kwargs: `VLineTrack` constructor kwargs; shorthand composition routes
                aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("vline", x_position=x_position, **kwargs)

    def _apply_theme_to_track(self, track: Track) -> Track:
        if self.theme is None:
            return track

        theme = self.theme
        for field_name in ("color", "alpha", "linewidth", "font_size", "cmap"):
            if not track.has_aesthetic(field_name):
                continue
            theme_value = getattr(theme, field_name)
            if theme_value is None:
                continue
            aesthetics_model = track.aesthetics_model()
            if aesthetics_model is None:
                continue
            model_field = aesthetics_model.model_fields.get(field_name)
            if model_field is None:
                continue
            if getattr(track, field_name) == model_field.default:
                setattr(track, field_name, theme_value)

        if hasattr(track, "label") and track.label is not None:
            default_label = track.label.__class__()
            for field_name, model_field in track.label.__class__.model_fields.items():
                current_value = getattr(track.label, field_name)
                theme_value = getattr(theme.label, field_name)
                default_value = getattr(default_label, field_name)
                if current_value == default_value and theme_value != default_value:
                    setattr(track.label, field_name, theme_value)

        return track

    @classmethod
    def available_track_aliases(cls) -> dict[str, str]:
        """Return available alias -> TrackClass mappings."""
        return {
            alias: track_cls.__name__
            for alias, track_cls in cls._alias_map().items()
        }

    @classmethod
    def track_options(cls, track: str | type[Track]) -> dict[str, dict]:
        """Return programmatic option metadata for a track alias or class.

        Examples:
            GenomicFigure.track_options("bigwig")
            GenomicFigure.track_options(BigWigTrack)
        """
        if isinstance(track, str):
            key = track.lower()
            alias_map = cls._alias_map()
            if key not in alias_map:
                raise ValueError(
                    f"Unknown track alias: {track}. Available: {list(alias_map.keys())}"
                )
            track_cls = alias_map[key]
        else:
            track_cls = track
        return list_options(track_cls)

    @classmethod
    def track_options_markdown(cls, track: str | type[Track]) -> str:
        """Return markdown-formatted option tables for a track alias or class."""
        if isinstance(track, str):
            key = track.lower()
            alias_map = cls._alias_map()
            if key not in alias_map:
                raise ValueError(
                    f"Unknown track alias: {track}. Available: {list(alias_map.keys())}"
                )
            track_cls = alias_map[key]
        else:
            track_cls = track
        return track_cls.options_markdown()

    def autoscale(self, enable: bool = True) -> Self:
        """Enable or disable automatic y-axis autoscaling across tracks."""
        self._autoscale = enable
        return self

    def autocolor(self, palette: str = "tab10") -> Self:
        """Apply a matplotlib palette across tracks that expose a `color` field."""
        self._autocolor_palette = palette
        import matplotlib.colors as mcolors

        cmap = plt.get_cmap(palette)
        color_index = 0
        for track in self.tracks:
            if self._should_autocolor_track(track):
                track.color = mcolors.to_hex(cmap(color_index % cmap.N))
                color_index += 1
        return self

    def highlight(self, region: str | GenomicRegion) -> Self:
        """Register a genomic region to draw as a background highlight."""
        self._highlight_regions.append(GenomicRegion.into(region))
        return self

    def highlight_style(self, color: str | None = None, alpha: float | None = None) -> Self:
        """Set default color/alpha used by highlighted regions."""
        if color is not None:
            self.highlight_color = color
        if alpha is not None:
            self.highlight_alpha = alpha
        return self

    def _create_track_from_alias(self, alias: str, **kwargs: Any) -> Track:
        alias_map = self._alias_map()
        key = alias.lower()
        if key not in alias_map:
            raise ValueError(
                f"Unknown track alias: {alias}. Available: {list(alias_map.keys())}"
            )
        track_cls = alias_map[key]

        track_fields = set(track_cls.model_fields.keys())
        aesthetics_model = track_cls.aesthetics_model()
        aesthetics_fields = (
            set(aesthetics_model.model_fields.keys()) if aesthetics_model is not None else set()
        )
        label_fields = set(LabelConfig.model_fields.keys())

        track_kwargs: dict[str, Any] = {}
        aesthetics_overrides: dict[str, Any] = {}
        label_overrides: dict[str, Any] = {}

        for field_name, field_value in kwargs.items():
            if field_name in {"aesthetics", "label"} or field_name in track_fields:
                track_kwargs[field_name] = field_value
            elif field_name in aesthetics_fields:
                aesthetics_overrides[field_name] = field_value
            elif field_name in label_fields:
                label_overrides[field_name] = field_value
            else:
                track_kwargs[field_name] = field_value

        if aesthetics_model is not None and aesthetics_overrides:
            base_aesthetics = track_kwargs.get("aesthetics")
            if isinstance(base_aesthetics, BaseModel):
                payload = base_aesthetics.model_dump()
            elif isinstance(base_aesthetics, dict):
                payload = dict(base_aesthetics)
            elif base_aesthetics is None:
                payload = {}
            else:
                payload = {}
            payload.update(aesthetics_overrides)
            track_kwargs["aesthetics"] = aesthetics_model(**payload)

        if label_overrides:
            base_label = track_kwargs.get("label")
            if isinstance(base_label, BaseModel):
                label_payload = base_label.model_dump()
            elif isinstance(base_label, dict):
                label_payload = dict(base_label)
            elif base_label is None:
                label_payload = {}
            else:
                label_payload = {}
            label_payload.update(label_overrides)
            track_kwargs["label"] = LabelConfig(**label_payload)

        return track_cls(**track_kwargs)

    @staticmethod
    def _alias_map() -> dict[str, type[Track]]:
        return {
            "scalebar": ScaleBar,
            "scale": ScaleBar,
            "genes": Genes,
            "spacer": Spacer,
            "bigwig": BigWigTrack,
            "bed": BedTrack,
            "axis": GenomicAxis,
            "highlight": HighlightsFromFile,
            "bigwig_overlay": OverlayTrack,
            "overlay": OverlayTrack,
            "narrowpeak": NarrowPeakTrack,
            "links": LinksTrack,
            "hline": HLineTrack,
            "vline": VLineTrack,
            "cooler": CoolerTrack,
            "capcruncher": CapcruncherTrack,
            "cooler_average": CoolerAverage,
            "bigwig_collection": BigWigCollection,
            "bigwig_diff": BigWigDiff,
        }

    @staticmethod
    def _apply_extend(gr: GenomicRegion, extend: float | int | None) -> GenomicRegion:
        if not extend:
            return gr

        if isinstance(extend, float):
            if extend < 0:
                raise ValueError("extend as float must be >= 0")
            bp = int(gr.length * extend)
        else:
            bp = int(extend)

        if bp < 0:
            raise ValueError("extend must be >= 0")
        return gr.extend(upstream=bp, downstream=bp)

    def _apply_autoscale_groups(
        self, axes: list[matplotlib.axes.Axes], main_tracks: list[Track]
    ) -> None:
        grouped: dict[str, list[int]] = {}
        for index, track in enumerate(main_tracks):
            if track.autoscale_group:
                grouped.setdefault(track.autoscale_group, []).append(index)

        for indices in grouped.values():
            mins = []
            maxs = []
            for index in indices:
                y_limits = axes[index].get_ylim()
                mins.append(y_limits[0])
                maxs.append(y_limits[1])
            group_min = min(mins)
            group_max = max(maxs)
            if group_min == group_max:
                group_max = group_min + 1
            for index in indices:
                axes[index].set_ylim(group_min, group_max)

    def _font_family(self) -> str:
        if self.theme is not None:
            if self.theme.font_family:
                return self.theme.font_family
            return self.theme.label.title_font
        return "DejaVu Sans"

    def plot(
        self,
        region: str | GenomicRegion,
        show: bool = False,
        extend: float | int | None = None,
        **kwargs: Any,
    ) -> matplotlib.figure.Figure | None:
        """Render all tracks for a single genomic region.

        Args:
            region: Region string (`chr:start-end`) or `GenomicRegion`.
            show: Whether to call `plt.show()` before returning.
            extend: Optional symmetric extension (fraction or base pairs).
            **kwargs: Forwarded to `matplotlib.figure.Figure` creation.

        Returns:
            The rendered matplotlib figure, or `None` if no tracks are present.
        """
        gr = GenomicRegion.into(region)
        gr = self._apply_extend(gr, extend)

        if not self.tracks:
            logger.warning("No tracks to plot")
            return None

        main_tracks = [track for track in self.tracks if track.height > 0]
        global_tracks = [track for track in self.tracks if track.height == 0]
        heights = [track.height for track in main_tracks]
        total_height = max(1.0, sum(heights) * self.track_height)

        fig = plt.figure(figsize=(self.width, total_height), **kwargs)
        fig.patch.set_facecolor("white")

        from matplotlib.gridspec import GridSpec

        gs = GridSpec(len(main_tracks), 1, height_ratios=heights, hspace=0.1)
        axes = [fig.add_subplot(gs[index]) for index in range(len(main_tracks))]

        def create_overlay_ax(zorder: int, label: str):
            ax = fig.add_axes(axes[0].get_position(), label=label, zorder=zorder)
            pos_top = axes[0].get_position()
            pos_bottom = axes[-1].get_position()
            ax.set_position(
                [
                    pos_top.x0,
                    pos_bottom.y0,
                    pos_top.width,
                    pos_top.y1 - pos_bottom.y0,
                ]
            )
            ax.patch.set_alpha(0)
            from .tracks.utils import clean_axis

            clean_axis(ax)
            ax.set_xlim(gr.start, gr.end)
            return ax

        bg_ax = create_overlay_ax(zorder=-1, label="background_overlay")
        fg_ax = create_overlay_ax(zorder=10, label="foreground_overlay")

        if self._autoscale:
            from .tracks.scaling import Autoscaler

            Autoscaler(tracks=self.tracks, gr=gr).apply()

        for ax, track in zip(axes, main_tracks):
            ax.patch.set_alpha(0)
            try:
                track.plot(ax, gr)
            except Exception as exc:
                logger.error(f"Error plotting track {track.__class__.__name__}: {exc}")
                ax.text(0.5, 0.5, f"Error: {exc}", ha="center", va="center", transform=ax.transAxes)

        self._apply_autoscale_groups(axes, main_tracks)

        for highlight_gr in self._highlight_regions:
            if highlight_gr.chromosome == gr.chromosome:
                start = max(highlight_gr.start, gr.start)
                end = min(highlight_gr.end, gr.end)
                if start < end:
                    bg_ax.axvspan(
                        start,
                        end,
                        color=self.highlight_color,
                        alpha=self.highlight_alpha,
                        zorder=-1,
                    )

        for track in global_tracks:
            try:
                if isinstance(track, HighlightsFromFile):
                    track.plot(bg_ax, gr)
                elif isinstance(track, VLineTrack):
                    track.plot(fg_ax, gr)
                elif hasattr(track, "plot_on_axes"):
                    track.plot_on_axes(gr, axes)
            except Exception as exc:
                logger.error(f"Error plotting global track {track.__class__.__name__}: {exc}")

        font_family = self._font_family()
        for axis in fig.axes:
            for text_artist in axis.texts:
                text_artist.set_fontfamily(font_family)

        fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)

        if show:
            plt.show()
        else:
            # Detach from pyplot state to avoid duplicate auto-rendering in notebooks
            # when the returned Figure object is also displayed as cell output.
            plt.close(fig)
        return fig

    def _resolve_gene_track(self) -> Genes:
        for track in self.tracks:
            if isinstance(track, Genes):
                return track
        return Genes(genome="hg38")

    def plot_gene(self, gene: str, extend: float = 0.5, **kwargs: Any) -> matplotlib.figure.Figure | None:
        """Resolve a gene symbol to a region and plot it.

        Args:
            gene: Gene symbol to look up in configured gene annotations.
            extend: Fractional extension around gene bounds.
            **kwargs: Forwarded to `plot`.

        Returns:
            The rendered figure or `None` when plotting is skipped.
        """
        genes_track = self._resolve_gene_track()
        has_genes_track = any(isinstance(track, Genes) for track in self.tracks)

        if genes_track.genome:
            bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
            with open(bed_prefix / "genes.json") as handle:
                mapping = json.load(handle)
            gene_file = bed_prefix / mapping[genes_track.genome]
            genes_df = pd.read_csv(gene_file, sep="\t", header=None)
        elif genes_track.data:
            genes_df = pd.read_csv(str(genes_track.data), sep="\t", header=None)
        else:
            raise ValueError("No gene annotation source available for plot_gene")

        genes_df.columns = [
            "chrom",
            "start",
            "end",
            "name",
            *[f"field_{index}" for index in range(max(0, genes_df.shape[1] - 4))],
        ]
        match = genes_df.loc[genes_df["name"].astype(str).str.upper() == gene.upper()]
        if match.empty:
            raise ValueError(f"Gene {gene} not found in annotation source")

        row = match.iloc[0]
        region = GenomicRegion(chromosome=row["chrom"], start=int(row["start"]), end=int(row["end"]))

        if has_genes_track:
            return self.plot(region, extend=extend, **kwargs)

        # If no genes track is present, include one for this plotting call so
        # the gene context is visible by default.
        self.tracks.append(genes_track)
        try:
            return self.plot(region, extend=extend, **kwargs)
        finally:
            self.tracks.pop()

    def plot_regions(
        self, regions: list[str] | str, ncols: int = 1, **kwargs: Any
    ) -> list[matplotlib.figure.Figure]:
        """Plot one or many regions, optionally composing a multi-column grid.

        Args:
            regions: A region string, list of region strings, or BED-like path.
            ncols: Number of columns for grid composition when >1 region.
            **kwargs: Forwarded to `plot`.

        Returns:
            A list of matplotlib figures.
        """
        if ncols < 1:
            raise ValueError("ncols must be >= 1")

        show = kwargs.pop("show", True)
        region_strings: list[str]
        path_candidate = Path(regions) if isinstance(regions, str) else None

        if isinstance(regions, str) and path_candidate.exists():
            bed_df = pd.read_csv(path_candidate, sep="\t", header=None, comment="#")
            region_strings = [
                f"{row[0]}:{int(row[1])}-{int(row[2])}" for _, row in bed_df.iterrows()
            ]
        elif isinstance(regions, str):
            region_strings = [regions]
        else:
            region_strings = regions

        figs: list[matplotlib.figure.Figure] = []
        for region in region_strings:
            fig = self.plot(region, show=False, **kwargs)
            if fig is not None:
                figs.append(fig)

        if not figs:
            return []

        if ncols > 1 and len(figs) > 1:
            nrows = math.ceil(len(figs) / ncols)
            grid_fig, grid_axes = plt.subplots(
                nrows,
                ncols,
                figsize=(self.width * ncols, self.track_height * max(2, nrows * 2)),
            )
            axes_array = np.atleast_1d(grid_axes).ravel()

            for index, source_fig in enumerate(figs):
                source_fig.canvas.draw()
                image = np.asarray(source_fig.canvas.buffer_rgba())
                axes_array[index].imshow(image)
                axes_array[index].set_axis_off()
                axes_array[index].set_title(region_strings[index], fontsize=9)

            for index in range(len(figs), len(axes_array)):
                axes_array[index].set_axis_off()

            for source_fig in figs:
                plt.close(source_fig)

            figs = [grid_fig]

        if show:
            for fig in figs:
                fig.show()
        return figs

    @staticmethod
    def _track_registry() -> dict[str, type[Track]]:
        return {
            cls.__name__: cls
            for cls in [
                BedTrack,
                BigWigTrack,
                BigWigCollection,
                BigWigDiff,
                BigwigOverlay,
                OverlayTrack,
                Genes,
                GenomicAxis,
                HighlightsFromFile,
                HLineTrack,
                LinksTrack,
                NarrowPeakTrack,
                ScaleBar,
                Spacer,
                VLineTrack,
                CoolerTrack,
                CapcruncherTrack,
                CoolerAverage,
            ]
        }

    def to_toml(self, path: str) -> None:
        """Serialize figure and tracks to a TOML file."""
        def _prune_none(value):
            if isinstance(value, dict):
                return {
                    key: _prune_none(item)
                    for key, item in value.items()
                    if item is not None
                }
            if isinstance(value, list):
                return [_prune_none(item) for item in value if item is not None]
            return value

        tracks_by_type: dict[str, list[dict]] = {}
        for track in self.tracks:
            tracks_by_type.setdefault(track.__class__.__name__, []).append(
                _prune_none(track.model_dump())
            )

        payload = {
            "figure": {"width": self.width, "track_height": self.track_height},
            "tracks": tracks_by_type,
        }
        payload = _prune_none(payload)

        try:
            import tomli_w
        except ImportError as exc:
            raise ImportError("to_toml requires optional dependency 'tomli-w'") from exc

        with open(path, "wb") as handle:
            handle.write(tomli_w.dumps(payload).encode("utf-8"))

    @classmethod
    def from_toml(cls, path: str) -> "GenomicFigure":
        """Load a `GenomicFigure` definition from a TOML file."""
        try:
            import tomllib
        except ImportError:
            import tomli as tomllib  # type: ignore[no-redef]

        with open(path, "rb") as handle:
            payload = tomllib.load(handle)

        figure_data = payload.get("figure", {})
        fig = cls(
            width=figure_data.get("width", 12),
            track_height=figure_data.get("track_height", 2.0),
        )

        registry = cls._track_registry()
        tracks_payload = payload.get("tracks", [])

        if isinstance(tracks_payload, list):
            for track_spec in tracks_payload:
                type_name = track_spec["type"]
                params = track_spec.get("params", {})
                if type_name not in registry:
                    raise ValueError(f"Unknown track type in TOML: {type_name}")
                fig.add_track(registry[type_name](**params))
            return fig

        if isinstance(tracks_payload, dict):
            for type_name, raw_params in tracks_payload.items():
                if type_name not in registry:
                    raise ValueError(f"Unknown track type in TOML: {type_name}")

                params_list = raw_params if isinstance(raw_params, list) else [raw_params]
                for params in params_list:
                    fig.add_track(registry[type_name](**(params or {})))
            return fig

        raise ValueError("Invalid TOML format: 'tracks' must be a list or table")

    def save(
        self,
        path: str | Path,
        region: str | GenomicRegion,
        dpi: int = 300,
        **kwargs: Any,
    ) -> None:
        """Render one region and write it to disk."""
        fig = self.plot(region, show=False)
        if fig:
            fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
            logger.info(f"Saved figure to {path}")
            plt.close(fig)

    def __repr__(self) -> str:
        parts = [f"tracks={len(self.tracks)}"]
        if self.width != 12:
            parts.append(f"width={self.width}")
        if self.theme is None:
            parts.append("theme=none")
        elif self.theme == Theme.minimal():
            parts.append("theme=minimal")
        elif self.theme == Theme.publication():
            parts.append("theme=publication")
        elif self.theme == Theme.default():
            parts.append("theme=default")
        else:
            parts.append("theme=custom")
        parts.append(f"autoscale={self._autoscale}")
        parts.append(f"highlights={len(self._highlight_regions)}")
        return f"GenomicFigure({', '.join(parts)})"

    def _repr_html_(self) -> str:
        import html

        if not self.tracks:
            return (
                "<div><em>No tracks added yet. Use .bigwig(), .genes(), .axis(), etc.</em></div>"
            )

        rows = []
        for index, track in enumerate(self.tracks):
            title = getattr(track, "title", "") or ""
            data_value = getattr(track, "data", None)
            data_display = "" if data_value is None else str(data_value)
            rows.append(
                "<tr>"
                f"<td>{index}</td>"
                f"<td>{html.escape(track.__class__.__name__)}</td>"
                f"<td>{html.escape(str(title))}</td>"
                f"<td>{html.escape(data_display)}</td>"
                f"<td>{html.escape(str(getattr(track, 'height', '')))}</td>"
                "</tr>"
            )

        table = (
            "<table>"
            "<thead><tr><th>#</th><th>Type</th><th>Title</th><th>Data</th><th>Height</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody>"
            "</table>"
        )
        return table


def _format_method_option_docs(alias: str) -> str:
    options = GenomicFigure.track_options(alias)
    lines = [
        "Auto-generated options (authoritative):",
        "",
        "Shorthand composition:",
        "- Pass track fields directly in kwargs.",
        "- Pass aesthetics fields directly in kwargs (auto-packed into `aesthetics`).",
        "- Pass label fields directly in kwargs (auto-packed into `label`).",
        "",
    ]

    for section in ("track", "aesthetics", "label"):
        section_data = options.get(section, {})
        lines.append(f"{section.title()} fields:")
        if not section_data:
            lines.append("- (none)")
            lines.append("")
            continue

        for field_name, meta in sorted(section_data.items()):
            description = meta.get("description") or "No description provided."
            choices = meta.get("choices") or []
            choices_text = f", choices={choices}" if choices else ""
            lines.append(
                f"- {field_name} ({meta['type']}, default={meta['default']}{choices_text}): {description}"
            )
        lines.append("")

    return "\n".join(lines).rstrip()


def _inject_figure_method_option_docs() -> None:
    method_alias_map = {
        "bigwig": "bigwig",
        "genes": "genes",
        "axis": "axis",
        "scalebar": "scalebar",
        "spacer": "spacer",
        "bed": "bed",
        "cooler": "cooler",
        "bigwig_collection": "bigwig_collection",
        "bigwig_diff": "bigwig_diff",
        "bigwig_overlay": "bigwig_overlay",
        "overlay": "overlay",
        "narrowpeak": "narrowpeak",
        "links": "links",
        "highlights": "highlight",
        "hline": "hline",
        "vline": "vline",
    }

    marker = "Auto-generated options (authoritative):"
    for method_name, alias in method_alias_map.items():
        method = getattr(GenomicFigure, method_name, None)
        if method is None:
            continue

        base_doc = (method.__doc__ or "").rstrip()
        if marker in base_doc:
            base_doc = base_doc.split(marker, 1)[0].rstrip()

        generated = _format_method_option_docs(alias)
        method.__doc__ = f"{base_doc}\n\n{generated}".strip()


_inject_figure_method_option_docs()
