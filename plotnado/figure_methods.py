"""GenomicFigureMethods — fluent track-builder mixin for GenomicFigure.

This module contains all the track-builder methods and their @overload type
stubs. It is separated from figure.py to keep the core class (plot rendering,
layout, theme application) readable.

Do not add rendering logic here. Every method must delegate to
``self.add_track(alias, ...)``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, Self, Unpack, overload

import pandas as pd

from ._kwargs import (
    AxisKwargs,
    BedKwargs,
    BigwigCollectionKwargs,
    BigwigDiffKwargs,
    BigwigKwargs,
    BigwigOverlayKwargs,
    CapcruncherKwargs,
    CoolerKwargs,
    CoolerAverageKwargs,
    GenesKwargs,
    HighlightsKwargs,
    HlineKwargs,
    LinksKwargs,
    NarrowpeakKwargs,
    OverlayKwargs,
    QuantnadoCoverageKwargs,
    QuantnadoMethylationKwargs,
    QuantnadoStrandedCoverageKwargs,
    QuantnadoVariantKwargs,
    ScalebarKwargs,
    SpacerKwargs,
    VlineKwargs,
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


class GenomicFigureMethods:
    """Fluent track-builder methods for adding tracks to a GenomicFigure.

    Mixed into ``GenomicFigure`` so that the core class stays focused on
    layout and rendering logic.
    """

    # BEGIN AUTO-GENERATED OVERLOAD: bigwig
    @overload
    def bigwig(
        self,
        data: Any,
        /,
        *,
        autoscale_group: str | None = ...,
        color_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        y_max: float | None = ...,
        y_min: float | None = ...,
        alpha: float = ...,
        baseline_alpha: float = ...,
        baseline_color: str = ...,
        baseline_linewidth: float = ...,
        color: str = ...,
        fill: bool = ...,
        linewidth: float = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        scatter_point_size: float = ...,
        show_baseline: bool = ...,
        smoothing_center: bool = ...,
        smoothing_method: Literal['mean', 'median'] = ...,
        smoothing_window: int = ...,
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
        color_group: str | None = ...,
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
        linewidth: float = ...,
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
        color_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        show_chromosome: bool = ...,
        title: str = ...,
        alpha: float = ...,
        axis_linewidth: float = ...,
        chromosome_fontweight: FontWeight = ...,
        color: str = ...,
        font_size: int = ...,
        linewidth: float = ...,
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
        color_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str = ...,
        alpha: float = ...,
        bar_linewidth: float = ...,
        color: str = ...,
        font_size: int = ...,
        label_offset: float = ...,
        linewidth: float = ...,
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
        color_group: str | None = ...,
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
        color_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        color: str = ...,
        display: DisplayMode = ...,
        draw_edges: bool = ...,
        edge_color: str = ...,
        font_size: int = ...,
        interval_height: float = ...,
        label_field: str = ...,
        linewidth: float = ...,
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
        color_group: str | None = ...,
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

    # BEGIN AUTO-GENERATED OVERLOAD: capcruncher
    @overload
    def capcruncher(
        self,
        file: str,
        /,
        *,
        autoscale_group: str | None = ...,
        balance: bool = ...,
        color_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        normalisation: str | None = ...,
        resolution: int | None = ...,
        title: str | None = ...,
        transform: CoolerTransform = ...,
        viewpoint: str | None = ...,
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
    # END AUTO-GENERATED OVERLOAD: capcruncher
    def capcruncher(self, file: str, /, **kwargs: Unpack[CapcruncherKwargs]) -> Self:
        """Add a capture-centric cooler matrix track.

        Args:
            file: Cooler file path.
            **kwargs: `CapcruncherTrack` constructor kwargs; shorthand composition
                routes aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("capcruncher", file=file, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: cooler_average
    @overload
    def cooler_average(
        self,
        files: list[str],
        /,
        *,
        autoscale_group: str | None = ...,
        balance: bool = ...,
        color_group: str | None = ...,
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
    # END AUTO-GENERATED OVERLOAD: cooler_average
    def cooler_average(
        self, files: list[str], /, **kwargs: Unpack[CoolerAverageKwargs]
    ) -> Self:
        """Add an averaged matrix track from multiple cooler files.

        Args:
            files: Cooler file paths to average.
            **kwargs: `CoolerAverage` constructor kwargs; shorthand composition
                routes aesthetics and label fields automatically.

        Returns:
            Self for method chaining.
        """
        return self.add_track("cooler_average", files=files, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: bigwig_collection
    @overload
    def bigwig_collection(
        self,
        files: list[str],
        /,
        *,
        autoscale_group: str | None = ...,
        color_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        colors: list[str] | None = ...,
        labels: list[str] | None = ...,
        linewidth: float = ...,
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
        color_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        method: BigWigDiffMethod = ...,
        title: str | None = ...,
        alpha: float = ...,
        bar_alpha: float = ...,
        color: str = ...,
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
        color_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        colors: list[str] | None = ...,
        linewidth: float = ...,
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
        color_group: str | None = ...,
        data: Any | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        colors: list[str] | None = ...,
        linewidth: float = ...,
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
        color_group: str | None = ...,
        height: float = ...,
        title: str | None = ...,
        alpha: float = ...,
        cmap: str = ...,
        color: str = ...,
        color_by: NarrowPeakColorBy | None = ...,
        display: DisplayMode = ...,
        draw_edges: bool = ...,
        edge_color: str = ...,
        font_size: int = ...,
        interval_height: float = ...,
        label_field: str = ...,
        linewidth: float = ...,
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
        color_group: str | None = ...,
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
        color_group: str | None = ...,
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
        color_group: str | None = ...,
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
        color_group: str | None = ...,
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

    # BEGIN AUTO-GENERATED OVERLOAD: quantnado_coverage
    @overload
    def quantnado_coverage(
        self,
        sample: str,
        /,
        *,
        autoscale_group: str | None = ...,
        color_group: str | None = ...,
        coverage_data: Any | None = ...,
        data: Any | None = ...,
        dataset_path: str | None = ...,
        height: float = ...,
        quantnado: Any | None = ...,
        title: str | None = ...,
        alpha: float = ...,
        baseline_alpha: float = ...,
        baseline_color: str = ...,
        baseline_linewidth: float = ...,
        color: str = ...,
        fill: bool = ...,
        linewidth: float = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        show_baseline: bool = ...,
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
    # END AUTO-GENERATED OVERLOAD: quantnado_coverage
    def quantnado_coverage(
        self, sample: str, /, **kwargs: Unpack[QuantnadoCoverageKwargs]
    ) -> Self:
        """Add a QuantNado coverage track for one sample.

        Args:
            sample: QuantNado sample name to render.
            **kwargs: `QuantNadoCoverageTrack` constructor kwargs.

        Returns:
            Self for method chaining.
        """
        return self.add_track("quantnado_coverage", sample=sample, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: quantnado_stranded_coverage
    @overload
    def quantnado_stranded_coverage(
        self,
        sample: str,
        /,
        *,
        autoscale_group: str | None = ...,
        color_group: str | None = ...,
        coverage_fwd_data: Any | None = ...,
        coverage_rev_data: Any | None = ...,
        data: Any | None = ...,
        dataset_path: str | None = ...,
        height: float = ...,
        quantnado: Any | None = ...,
        title: str | None = ...,
        alpha: float = ...,
        baseline_alpha: float = ...,
        baseline_color: str = ...,
        baseline_linewidth: float = ...,
        color: str = ...,
        fill: bool = ...,
        linewidth: float = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        reverse_color: str | None = ...,
        show_baseline: bool = ...,
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
    # END AUTO-GENERATED OVERLOAD: quantnado_stranded_coverage
    def quantnado_stranded_coverage(
        self, sample: str, /, **kwargs: Unpack[QuantnadoStrandedCoverageKwargs]
    ) -> Self:
        """Add a QuantNado stranded coverage track for one sample.

        Args:
            sample: QuantNado sample name to render.
            **kwargs: `QuantNadoStrandedCoverageTrack` constructor kwargs.

        Returns:
            Self for method chaining.
        """
        return self.add_track("quantnado_stranded_coverage", sample=sample, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: quantnado_methylation
    @overload
    def quantnado_methylation(
        self,
        sample: str,
        /,
        *,
        autoscale_group: str | None = ...,
        color_group: str | None = ...,
        data: Any | None = ...,
        dataset_path: str | None = ...,
        height: float = ...,
        methylation_data: Any | None = ...,
        methylation_variable: str = ...,
        quantnado: Any | None = ...,
        title: str | None = ...,
        alpha: float = ...,
        color: str = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        point_size: float = ...,
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
    # END AUTO-GENERATED OVERLOAD: quantnado_methylation
    def quantnado_methylation(
        self, sample: str, /, **kwargs: Unpack[QuantnadoMethylationKwargs]
    ) -> Self:
        """Add a QuantNado methylation track for one sample.

        Args:
            sample: QuantNado sample name to render.
            **kwargs: `QuantNadoMethylationTrack` constructor kwargs.

        Returns:
            Self for method chaining.
        """
        return self.add_track("quantnado_methylation", sample=sample, **kwargs)

    # BEGIN AUTO-GENERATED OVERLOAD: quantnado_variant
    @overload
    def quantnado_variant(
        self,
        sample: str,
        /,
        *,
        allele_depth_alt_data: Any | None = ...,
        allele_depth_alt_variable: str = ...,
        allele_depth_ref_data: Any | None = ...,
        allele_depth_ref_variable: str = ...,
        autoscale_group: str | None = ...,
        color_group: str | None = ...,
        data: Any | None = ...,
        dataset_path: str | None = ...,
        fetch_genotype: bool = ...,
        genotype_data: Any | None = ...,
        genotype_variable: str = ...,
        height: float = ...,
        quantnado: Any | None = ...,
        title: str | None = ...,
        alpha: float = ...,
        baseline_alpha: float = ...,
        baseline_color: str = ...,
        baseline_linewidth: float = ...,
        het_color: str = ...,
        hom_alt_color: str = ...,
        linewidth: float = ...,
        marker_size: float = ...,
        max_value: float | None = ...,
        min_value: float | None = ...,
        show_baseline: bool = ...,
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
    # END AUTO-GENERATED OVERLOAD: quantnado_variant
    def quantnado_variant(
        self, sample: str, /, **kwargs: Unpack[QuantnadoVariantKwargs]
    ) -> Self:
        """Add a QuantNado variant allele-frequency track for one sample.

        Args:
            sample: QuantNado sample name to render.
            **kwargs: `QuantNadoVariantTrack` constructor kwargs.

        Returns:
            Self for method chaining.
        """
        return self.add_track("quantnado_variant", sample=sample, **kwargs)
