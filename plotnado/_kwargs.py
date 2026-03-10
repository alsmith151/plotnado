from __future__ import annotations

from pathlib import Path
from typing import Any, TypedDict

import pandas as pd

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
from .tracks.region import GenomicRegion
from .tracks.base import Track



class BigwigKwargs(TypedDict, total=False):
    title: str | None
    height: float
    autoscale_group: str | None
    y_min: float | None
    y_max: float | None
    style: PlotStyle
    color: str
    fill: bool
    alpha: float
    linewidth: float
    scatter_point_size: float
    show_baseline: bool
    baseline_color: str
    baseline_alpha: float
    baseline_linewidth: float
    smoothing_window: int
    smoothing_method: Literal['mean', 'median']
    smoothing_center: bool
    min_value: float | None
    max_value: float | None
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class GenesKwargs(TypedDict, total=False):
    title: str | None
    data: Path | str | DataFrame | None
    height: float
    autoscale_group: str | None
    row_scale: float
    small_relative: float
    gene_count: int
    style: PlotStyle
    color: str
    fill: bool
    alpha: float
    display: DisplayMode
    minimum_gene_length: int
    max_number_of_rows: int
    interval_height: float
    arrow_size: float
    exon_linewidth: float
    exon_edge_color: str
    exon_color: str
    intron_linewidth: float
    intron_color: str
    chevron_height_ratio: float
    chevron_vertical_offset_ratio: float
    chevron_width_fraction: float
    chevron_min_width_bp: float
    chevron_margin_bp: float
    chevron_target_spacing_bp: float
    chevron_max_count: int
    gene_label_font_size: int
    gene_label_style: GeneLabelStyle
    show_labels: bool
    label_overlap_strategy: GeneLabelOverlapStrategy
    label_max_chars: int
    label_offset_fraction: float
    label_stagger_offset: float
    label_vertical_offset: float
    label_min_overlap_bp: int
    label_min_overlap_fraction: float
    label_connectors: bool
    label_connector_linewidth: float
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class AxisKwargs(TypedDict, total=False):
    title: str
    data: Any | None
    height: float
    autoscale_group: str | None
    show_chromosome: bool
    color: str
    font_size: int
    num_ticks: int
    use_human_readable_labels: bool
    tick_height: float
    axis_linewidth: float
    tick_color: str
    tick_linewidth: float
    chromosome_fontweight: FontWeight
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class ScalebarKwargs(TypedDict, total=False):
    title: str
    data: Any | None
    height: float
    autoscale_group: str | None
    style: PlotStyle
    color: str
    position: Position
    scale_distance: float | None
    font_size: int
    bar_linewidth: float
    tick_linewidth: float
    tick_height: float
    label_offset: float
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class SpacerKwargs(TypedDict, total=False):
    title: str
    data: Any | None
    autoscale_group: str | None
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class BedKwargs(TypedDict, total=False):
    title: str | None
    height: float
    autoscale_group: str | None
    color: str
    edge_color: str
    alpha: float
    interval_height: float
    display: DisplayMode
    max_rows: int
    show_labels: bool
    label_field: str
    font_size: int
    rect_linewidth: float
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class CoolerKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    resolution: int | None
    balance: bool
    transform: CoolerTransform
    cmap: str
    min_value: float | None
    max_value: float | None
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class BigwigCollectionKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    colors: list[str] | None
    labels: list[str] | None
    alpha: float
    style: CollectionStyle
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class BigwigDiffKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    method: BigWigDiffMethod
    positive_color: str
    negative_color: str
    linewidth: float
    bar_alpha: float
    zero_line_color: str
    zero_line_width: float
    zero_line_alpha: float
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class BigwigOverlayKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    colors: list[str] | None
    alpha: float
    show_labels: bool
    min_value: float | None
    max_value: float | None
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class OverlayKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    colors: list[str] | None
    alpha: float
    show_labels: bool
    min_value: float | None
    max_value: float | None
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class NarrowpeakKwargs(TypedDict, total=False):
    title: str | None
    height: float
    autoscale_group: str | None
    color: str
    edge_color: str
    alpha: float
    interval_height: float
    display: DisplayMode
    max_rows: int
    show_labels: bool
    label_field: str
    font_size: int
    rect_linewidth: float
    color_by: NarrowPeakColorBy | None
    cmap: str
    min_score: float | None
    max_score: float | None
    show_summit: bool
    summit_color: str
    summit_width: float
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class LinksKwargs(TypedDict, total=False):
    title: str | None
    height: float
    autoscale_group: str | None
    color: str
    edge_color: str | None
    alpha: float
    linewidth: float
    cmap: str
    max_height: float
    color_by_score: bool
    min_score: float | None
    max_score: float | None
    y_baseline: float
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class HighlightsKwargs(TypedDict, total=False):
    title: str | None
    height: float
    autoscale_group: str | None
    color: str
    alpha: float
    edge_color: str | None
    linewidth: float
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class HlineKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    color: str
    linestyle: str
    linewidth: float
    alpha: float
    zorder: int
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class VlineKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    color: str
    linestyle: str
    linewidth: float
    alpha: float
    zorder: int
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class CapcruncherKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    resolution: int | None
    balance: bool
    transform: CoolerTransform
    viewpoint: str | None
    normalisation: str | None
    cmap: str
    min_value: float | None
    max_value: float | None
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight

class CoolerAverageKwargs(TypedDict, total=False):
    title: str | None
    data: Any | None
    height: float
    autoscale_group: str | None
    resolution: int | None
    balance: bool
    transform: CoolerTransform
    cmap: str
    min_value: float | None
    max_value: float | None
    plot_title: bool
    plot_scale: bool
    label_on_track: bool
    data_range_style: DataRangeStyle
    label_box_enabled: bool
    label_box_alpha: float
    title_location: Position
    title_height: float
    title_size: int
    title_color: str
    title_font: str
    title_weight: FontWeight
    scale_location: Position
    scale_height: float
    scale_precision: int
    scale_size: int
    scale_color: str
    scale_font: str
    scale_weight: FontWeight
