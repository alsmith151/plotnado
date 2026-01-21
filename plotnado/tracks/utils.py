"""
General utilities for genomic tracks.
"""

import matplotlib.axes
import matplotlib.ticker


def clean_axis(ax: matplotlib.axes.Axes) -> matplotlib.axes.Axes:
    """
    Clean the axis by removing ticks and spines.

    Args:
        ax: The matplotlib axis object to clean
    """
    ax.set_xticks([])
    ax.set_yticks([])
    if hasattr(ax.spines, "items"):
        for spine in ax.spines.values():
            spine.set_visible(False)
    else:
        if isinstance(ax.spines, dict):
            for key in ["top", "right", "left", "bottom"]:
                if key in ax.spines:
                    ax.spines[key].set_visible(False)

    if hasattr(ax.xaxis, "set_major_locator"):
        ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
    if hasattr(ax.yaxis, "set_major_locator"):
        ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())

    return ax


def get_human_readable_number_of_bp(bp: int) -> str:
    """Convert integer into human readable basepair number."""
    if bp < 1000:
        return f"{bp: .0f} bp"
    elif (bp / 1e3) < 1000:
        return f"{bp / 1e3: .0f} kb"
    elif (bp / 1e6) < 1000:
        return f"{bp / 1e6: .0f} mb"
    return f"{bp / 1e9: .0f} gb"
