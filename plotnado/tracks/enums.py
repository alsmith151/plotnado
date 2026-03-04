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
    """Available track types for Figure.add_track()."""

    BIGWIG = "bigwig"
    GENES = "genes"
    SCALEBAR = "scalebar"
    SCALE = "scale"  # Alias for scalebar
    SPACER = "spacer"
    BED = "bed"
    GENOMIC_AXIS = "axis"
    HIGHLIGHT = "highlight"
