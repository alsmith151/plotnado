"""
Enums for plotnado track configuration.

Provides clear, type-safe options for track aesthetics and behavior.
"""

from enum import Enum


class PlotStyle(str, Enum):
    """Style for plotting signal tracks."""

    STD = "std"
    FILL = "fill"
    LINE = "line"
    SCATTER = "scatter"
    HEATMAP = "heatmap"
    FRAGMENT = "fragment"


class CollectionStyle(str, Enum):
    """Layout style for collections of signal tracks."""

    OVERLAY = "overlay"
    STACKED = "stacked"


class Position(str, Enum):
    """Position options for elements within tracks."""

    LEFT = "left"
    RIGHT = "right"
    CENTER = "center"


class Strand(str, Enum):
    """DNA strand orientation."""

    PLUS = "+"
    MINUS = "-"


class DisplayMode(str, Enum):
    """Display mode for interval tracks like genes."""

    COLLAPSED = "collapsed"
    EXPANDED = "expanded"


class FontWeight(str, Enum):
    """Font weight options."""

    NORMAL = "normal"
    BOLD = "bold"


class DataRangeStyle(str, Enum):
    """How y-range information is displayed on tracks."""

    TEXT = "text"
    COLORBAR = "colorbar"
    NONE = "none"


class GeneLabelStyle(str, Enum):
    """Font style options for gene labels."""

    NORMAL = "normal"
    ITALIC = "italic"
    OBLIQUE = "oblique"


class GeneLabelOverlapStrategy(str, Enum):
    """Collision handling strategies for gene labels."""

    AUTO = "auto"
    SMART = "smart"
    STAGGER = "stagger"
    SUPPRESS = "suppress"
    AUTO_EXPAND = "auto_expand"


class CoolerTransform(str, Enum):
    """Matrix transform options for cooler-like tracks."""

    LOG = "log"
    LOG2 = "log2"
    LOG10 = "log10"
    NONE = "none"


class BigWigDiffMethod(str, Enum):
    """Computation method for differential BigWig tracks."""

    SUBTRACT = "subtract"
    RATIO = "ratio"
    LOG2RATIO = "log2ratio"


class NarrowPeakColorBy(str, Enum):
    """Fields that can drive narrowPeak color mapping."""

    SCORE = "score"
    SIGNAL_VALUE = "signalValue"
    P_VALUE = "pValue"
    Q_VALUE = "qValue"


class ScalingMethod(str, Enum):
    """Supported scaling factor aggregation methods."""

    MAX = "max"
    MEAN = "mean"
    TOTAL = "total"


class TrackType(str, Enum):
    """All available track types, used in both the Python API and YAML templates.

    Values match the string aliases accepted by ``GenomicFigure.add_track()``
    and the ``type`` field in YAML template files.

    Example:
        >>> from plotnado.tracks.enums import TrackType
        >>> TrackType.BIGWIG
        <TrackType.BIGWIG: 'bigwig'>
        >>> TrackType("bigwig")
        <TrackType.BIGWIG: 'bigwig'>
    """

    # Signal tracks
    BIGWIG = "bigwig"
    BEDGRAPH = "bedgraph"          # Rendered by BigWigTrack
    BIGWIG_DIFF = "bigwig_diff"
    BIGWIG_COLLECTION = "bigwig_collection"
    BIGWIG_OVERLAY = "bigwig_overlay"

    # Interval tracks
    BED = "bed"
    ANNOTATION = "annotation"      # Semantic alias for BED
    NARROWPEAK = "narrowpeak"
    LINKS = "links"

    # Gene annotation
    GENE = "gene"
    GENES = "genes"                # Alias accepted by add_track()

    # Guide tracks
    SCALEBAR = "scalebar"
    SCALE = "scale"                # Alias for scalebar
    AXIS = "axis"
    SPACER = "spacer"

    # Overlay
    OVERLAY = "overlay"

    # Highlight / annotation lines
    HIGHLIGHT = "highlight"
    HLINE = "hline"
    VLINE = "vline"

    # Hi-C / contact matrices
    COOLER = "cooler"
    CAPCRUNCHER = "capcruncher"
    COOLER_AVERAGE = "cooler_average"

    # QuantNado multi-omics
    QUANTNADO_COVERAGE = "quantnado_coverage"
    QUANTNADO_STRANDED_COVERAGE = "quantnado_stranded_coverage"
    QUANTNADO_METHYLATION = "quantnado_methylation"
    QUANTNADO_VARIANT = "quantnado_variant"

    # Fallback
    UNKNOWN = "unknown"            # Falls back to BED rendering
