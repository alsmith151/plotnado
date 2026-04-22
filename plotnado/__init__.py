"""
Plotnado - Simple genomic track visualization

A lightweight Python package for creating genome browser-style plots
without heavy dependencies like CoolBox.
"""

from .figure import GenomicFigure
from .igv import parse_igv_session, IgvSession
from .theme import Theme
from .template import Template, TrackSpec, GuideSpec, GroupSpec
from .render import TemplateCompiler, RenderPlan, ResolvedTrack
from .tracks import (
    # Enums
    CollectionStyle,
    DisplayMode,
    FontWeight,
    PlotStyle,
    Position,
    Strand,
    TrackType,
    # Base classes
    GenomicRegion,
    LabelConfig,
    Track,
    TrackLabeller,
    list_options,
    # Track types
    BigWigTrack,
    BigwigAesthetics,
    Genes,
    GeneTrack,
    GenesAesthetics,
    ScaleBar,
    ScaleBarTrack,
    ScaleBarAesthetics,
    Spacer,
    SpacerTrack,
    BedTrack,
    BedAesthetics,
    GenomicAxis,
    AxisTrack,
    GenomicAxisAesthetics,
    HighlightsFromFile,
    HighlightsAesthetics,
    OverlayTrack,
    OverlayTrackAesthetics,
    BigwigOverlay,
    BigwigOverlayAesthetics,
    BigWigCollection,
    BigWigCollectionAesthetics,
    BigWigDiff,
    BigWigDiffAesthetics,
    CoolerTrack,
    CapcruncherTrack,
    CoolerAverage,
    NarrowPeakTrack,
    NarrowPeakAesthetics,
    LinksTrack,
    LinksAesthetics,
    HLineTrack,
    VLineTrack,
    AnnotationAesthetics,
    Autoscaler,
    QuantNadoCoverageTrack,
    QuantNadoCoverageAesthetics,
    QuantNadoStrandedCoverageTrack,
    QuantNadoStrandedCoverageAesthetics,
    QuantNadoMethylationTrack,
    QuantNadoMethylationAesthetics,
    QuantNadoVariantTrack,
    QuantNadoVariantAesthetics,
)

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

__all__ = [
    "GenomicFigure",
    "parse_igv_session",
    "IgvSession",
    "Theme",
    # Template system
    "Template",
    "TrackSpec",
    "GuideSpec",
    "GroupSpec",
    # Render pipeline
    "TemplateCompiler",
    "RenderPlan",
    "ResolvedTrack",
    # Enums
    "DisplayMode",
    "CollectionStyle",
    "FontWeight",
    "PlotStyle",
    "Position",
    "Strand",
    "TrackType",
    # Base
    "GenomicRegion",
    "LabelConfig",
    "Track",
    "TrackLabeller",
    "list_options",
    # Track types
    "BigWigTrack",
    "BigwigAesthetics",
    "Genes",
    "GeneTrack",
    "GenesAesthetics",
    "ScaleBar",
    "ScaleBarTrack",
    "ScaleBarAesthetics",
    "Spacer",
    "SpacerTrack",
    "BedTrack",
    "BedAesthetics",
    "GenomicAxis",
    "AxisTrack",
    "GenomicAxisAesthetics",
    "HighlightsFromFile",
    "HighlightsAesthetics",
    "OverlayTrack",
    "OverlayTrackAesthetics",
    "BigwigOverlay",
    "BigwigOverlayAesthetics",
    "BigWigCollection",
    "BigWigCollectionAesthetics",
    "BigWigDiff",
    "BigWigDiffAesthetics",
    "CoolerTrack",
    "CapcruncherTrack",
    "CoolerAverage",
    "NarrowPeakTrack",
    "NarrowPeakAesthetics",
    "LinksTrack",
    "LinksAesthetics",
    "HLineTrack",
    "VLineTrack",
    "AnnotationAesthetics",
    "Autoscaler",
    "QuantNadoCoverageTrack",
    "QuantNadoCoverageAesthetics",
    "QuantNadoStrandedCoverageTrack",
    "QuantNadoStrandedCoverageAesthetics",
    "QuantNadoMethylationTrack",
    "QuantNadoMethylationAesthetics",
    "QuantNadoVariantTrack",
    "QuantNadoVariantAesthetics",
]
