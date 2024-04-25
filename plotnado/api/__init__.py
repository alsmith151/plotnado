import coolbox.api as cb

from .tracks import (
    Autoscaler,
    BedMemory,
    BedSimple,
    BigwigFragment,
    BigwigFragmentCollection,
    BigwigFragmentCollectionOverlay,
    BigwigOverlay,
    GenomicAxis,
    HighlightsFromFile,
    MatrixCapcruncher,
    MatrixCapcruncherAverage,
    ScaleBar,
)
from .genes import Genes
from .track_wrapper import (
    BIGWIG_TRACKS,
    CUSTOM_TRACKS,
    FilelessTracks,
    MATRIX_TRACKS,
    TrackType,
    TrackWrapper,
)
from .figure import Figure


__all__ = [
    "Autoscaler",
    "BedMemory",
    "BedSimple",
    "BigwigFragment",
    "BigwigFragmentCollection",
    "BigwigFragmentCollectionOverlay",
    "BigwigOverlay",
    "Figure",
    "Genes",
    "GenomicAxis",
    "HighlightsFromFile",
    "MatrixCapcruncher",
    "MatrixCapcruncherAverage",
    "ScaleBar",
    "TrackWrapper",
]
