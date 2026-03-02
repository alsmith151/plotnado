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
from .region import GenomicRegion, GenomicRange
from .utils import (
    clean_axis,
    format_distance,
    get_human_readable_number_of_bp,
    format_genomic_value,
    parse_genomic_value,
)
from .base import (
    Track,
    TrackLabeller,
)
from .aesthetics import Aesthetics
from .schemas import BedgraphDataFrame
from .bigwig import BigWigTrack, BigwigAesthetics
from .genes import Genes, GenesAesthetics
from .scalebar import ScaleBar, ScaleBarAesthetics
from .spacer import Spacer
from .bed import BedTrack, BedAesthetics
from .axis import GenomicAxis, GenomicAxisAesthetics
from .highlight import HighlightsFromFile, HighlightsAesthetics
from .scaling import Autoscaler, Scaler
from .overlay import BigwigOverlay, BigwigOverlayAesthetics
from .bigwig_collection import BigWigCollection, BigWigCollectionAesthetics
from .bigwig_diff import BigWigDiff
from .cooler_track import CoolerTrack, CapcruncherTrack, CoolerAverage
from .peaks import NarrowPeakTrack, NarrowPeakAesthetics
from .links import LinksTrack, LinksAesthetics
from .annotations import HLineTrack, VLineTrack, AnnotationAesthetics

GeneTrack = Genes
AxisTrack = GenomicAxis
ScaleBarTrack = ScaleBar
SpacerTrack = Spacer

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
    "GenomicRange",
    "Track",
    "TrackLabeller",
    "clean_axis",
    "format_distance",
    "get_human_readable_number_of_bp",
    "format_genomic_value",
    "parse_genomic_value",
    # Data schemas
    "BedgraphDataFrame",
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
    "Autoscaler",
    "Scaler",
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
]
