"""
Plotnado - Simple genomic track visualization

A lightweight Python package for creating genome browser-style plots
without heavy dependencies like CoolBox.
"""

from .figure import Figure
from .tracks import (
    # Enums
    DisplayMode,
    FontWeight,
    PlotStyle,
    Position,
    Strand,
    TrackType,
    # Base classes
    GenomicRegion,
    Track,
    TrackLabeller,
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
    BigwigOverlay,
    BigwigOverlayAesthetics,
    BigWigCollection,
    BigWigCollectionAesthetics,
    BigWigDiff,
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
    # Enums
    "DisplayMode",
    "FontWeight",
    "PlotStyle",
    "Position",
    "Strand",
    "TrackType",
    # Base
    "GenomicRegion",
    "Track",
    "TrackLabeller",
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
    "BigwigOverlay",
    "BigwigOverlayAesthetics",
    "BigWigCollection",
    "BigWigCollectionAesthetics",
    "BigWigDiff",
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
