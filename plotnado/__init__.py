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
    GenesAesthetics,
    ScaleBar,
    ScaleBarAesthetics,
    Spacer,
    BedTrack,
    BedAesthetics,
    GenomicAxis,
    GenomicAxisAesthetics,
    HighlightsFromFile,
    HighlightsAesthetics,
    BigwigOverlay,
    BigwigOverlayAesthetics,
    Autoscaler,
    tabix_gtf,
)

__version__ = "0.3.0"

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
    "GenesAesthetics",
    "ScaleBar",
    "ScaleBarAesthetics",
    "Spacer",
    "BedTrack",
    "BedAesthetics",
    "GenomicAxis",
    "GenomicAxisAesthetics",
    "HighlightsFromFile",
    "HighlightsAesthetics",
    "BigwigOverlay",
    "BigwigOverlayAesthetics",
    "Autoscaler",
    "tabix_gtf",
]
