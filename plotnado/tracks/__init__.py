"""
Plotnado Tracks Module

Provides track classes for genomic visualization.
"""

from .enums import (
    CollectionStyle,
    DisplayMode,
    FontWeight,
    GeneLabelOverlapStrategy,
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
    LabelConfig,
    Track,
    TrackLabeller,
)
from .aesthetics import list_options
from .schemas import BedgraphDataFrame
from .bigwig import BigWigTrack, BigwigAesthetics
from .genes import Genes, GenesAesthetics
from .scalebar import ScaleBar, ScaleBarAesthetics
from .spacer import Spacer
from .bed import BedTrack, BedAesthetics
from .axis import GenomicAxis, GenomicAxisAesthetics
from .highlight import HighlightsFromFile, HighlightsAesthetics
from .scaling import Autoscaler, Scaler
from .overlay import OverlayTrack, OverlayTrackAesthetics, BigwigOverlay, BigwigOverlayAesthetics
from .bigwig_collection import BigWigCollection, BigWigCollectionAesthetics
from .bigwig_diff import BigWigDiff, BigWigDiffAesthetics
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
    "CollectionStyle",
    "FontWeight",
    "GeneLabelOverlapStrategy",
    "PlotStyle",
    "Position",
    "Strand",
    "TrackType",
    # Base classes
    "list_options",
    "GenomicRegion",
    "GenomicRange",
    "Track",
    "TrackLabeller",
    "LabelConfig",
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
]
