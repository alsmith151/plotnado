"""
Plotnado - Simple genomic track visualization

A lightweight Python package for creating genome browser-style plots
without heavy dependencies like CoolBox.
"""

from .figure import Figure
from .theme import Theme
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
)

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

__all__ = [
    "Figure",
    "Theme",
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
]
