"""
Plotnado Tracks Module

Provides track classes for genomic visualization.
"""

from .enums import (
    DisplayMode,
    FontWeight,
    PlotStyle,
    Position,
    Strand,
    TrackType,
)
from .region import GenomicRegion
from .utils import clean_axis, get_human_readable_number_of_bp
from .base import (
    Track,
    TrackLabeller,
)
from .aesthetics import Aesthetics
from .schemas import BedgraphDataFrame
from .bigwig import BigWigTrack, BigwigAesthetics
from .genes import Genes, GenesAesthetics, tabix_gtf
from .scalebar import ScaleBar, ScaleBarAesthetics
from .spacer import Spacer
from .bed import BedTrack, BedAesthetics
from .axis import GenomicAxis, GenomicAxisAesthetics
from .highlight import HighlightsFromFile, HighlightsAesthetics
from .scaling import Autoscaler, Scaler
from .overlay import BigwigOverlay, BigwigOverlayAesthetics

__all__ = [
    # Enums
    "DisplayMode",
    "FontWeight",
    "PlotStyle",
    "Position",
    "Strand",
    "TrackType",
    # Base classes
    "Aesthetics",
    "GenomicRegion",
    "Track",
    "TrackLabeller",
    "clean_axis",
    "get_human_readable_number_of_bp",
    # Data schemas
    "BedgraphDataFrame",
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
    "Autoscaler",
    "Scaler",
    "BigwigOverlay",
    "BigwigOverlayAesthetics",
    "tabix_gtf",
]
