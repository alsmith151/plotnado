"""
Enums for plotnado track configuration.

Provides clear, type-safe options for track aesthetics and behavior.
"""

from enum import Enum


class PlotStyle(str, Enum):
    """Style for plotting signal tracks."""

    FILL = "fill"
    LINE = "line"
    SCATTER = "scatter"
    HEATMAP = "heatmap"


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
