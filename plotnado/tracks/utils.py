"""
General utilities for genomic tracks.
"""

from functools import lru_cache
import matplotlib.axes
import matplotlib.ticker
import pandas as pd
import pybigtools


import re


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


def parse_genomic_value(value: str | int | float) -> int:
    """
    Parse a genomic value that might have suffixes (k, M, G) or commas.

    Examples:
        '1k' -> 1000
        '1.2M' -> 1,200,000
        '1,500,000' -> 1,500,000
        '20kb' -> 20,000
    """
    if isinstance(value, (int, float)):
        return int(value)

    # Remove commas
    s = str(value).replace(",", "").strip().lower()

    # Handle suffixes
    multipliers = {
        "k": 1_000,
        "kb": 1_000,
        "m": 1_000_000,
        "mb": 1_000_000,
        "g": 1_000_000_000,
        "gb": 1_000_000_000,
    }

    # Find suffix
    match = re.match(r"^([\d.]+)\s*([a-z]+)?$", s)
    if not match:
        raise ValueError(f"Invalid genomic value: {value}")

    num_str, suffix = match.groups()
    val = float(num_str)

    if suffix:
        if suffix in multipliers:
            val *= multipliers[suffix]
        else:
            raise ValueError(f"Unknown suffix: {suffix} in {value}")

    return int(val)


def format_genomic_value(pos: int, precision: int | None = None) -> str:
    """Format a genomic position for display (e.g. 1.2M, 450k).

    Args:
        pos: Genomic position in base pairs
        precision: Number of decimal places. If None, auto-detect minimal precision.
    """
    if pos >= 1_000_000:
        val = pos / 1_000_000
        if precision is not None:
            return f"{val:.{precision}f}M"
        if val % 1 == 0:
            return f"{val:.0f}M"
        return f"{val:.2f}M".rstrip("0").rstrip(".")
    elif pos >= 1_000:
        val = pos / 1_000
        if precision is not None:
            return f"{val:.{precision}f}k"
        if val % 1 == 0:
            return f"{val:.0f}k"
        return f"{val:.1f}k".rstrip("0").rstrip(".")
    return f"{pos:,}"


def format_distance(bp: int) -> str:
    """Convert integer into human readable basepair distance (e.g. 20 kb)."""
    if bp < 1000:
        return f"{bp: .0f} bp"
    elif (bp / 1e3) < 1000:
        val = bp / 1e3
        if val % 1 == 0:
            return f"{val:.0f} kb"
        return f"{val:.1f} kb".rstrip("0").rstrip(".")
    elif (bp / 1e6) < 1000:
        val = bp / 1e6
        if val % 1 == 0:
            return f"{val:.0f} Mb"
        return f"{val:.1f} Mb".rstrip("0").rstrip(".")
    return f"{bp / 1e9: .0f} Gb"


def get_human_readable_number_of_bp(bp: int) -> str:
    """Backward-compatible distance formatter used by historical tests."""
    if bp < 1000:
        return f" {int(bp)} bp"
    if (bp / 1e3) < 1000:
        return f" {int(round(bp / 1e3))} kb"
    if (bp / 1e6) < 1000:
        return f" {int(round(bp / 1e6))} mb"
    return f" {int(round(bp / 1e9))} gb"


def intervals_to_dataframe(records: list[tuple], chrom: str) -> pd.DataFrame:
    """Convert indexed tuple records from pybigtools into a DataFrame."""
    if not records:
        return pd.DataFrame(columns=["chrom", "start", "end"])

    max_width = max(len(record) for record in records)
    extra_fields = [f"field_{i}" for i in range(1, max_width - 1)]
    columns = ["start", "end", *extra_fields]
    padded = [tuple(record) + (None,) * (max_width - len(record)) for record in records]

    df = pd.DataFrame(padded, columns=columns)
    df.insert(0, "chrom", chrom)
    return df


def _normalise_interval_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    if "Chromosome" in df.columns:
        rename_map["Chromosome"] = "chrom"
    if "Start" in df.columns:
        rename_map["Start"] = "start"
    if "End" in df.columns:
        rename_map["End"] = "end"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def _slice_region(df: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    if df.empty:
        return df

    data = _normalise_interval_columns(df)
    if not {"chrom", "start", "end"}.issubset(data.columns):
        return pd.DataFrame(columns=["chrom", "start", "end"])

    return data.loc[
        (data["chrom"] == chrom) & (data["end"] > start) & (data["start"] < end)
    ].copy()


@lru_cache(maxsize=32)
def _read_bed_cached(path: str) -> pd.DataFrame:
    import pyranges1 as pr

    gr = pr.read_bed(path)
    if hasattr(gr, "df"):
        return gr.df.copy()
    if hasattr(gr, "as_df"):
        return gr.as_df().copy()
    return pd.DataFrame(gr)


@lru_cache(maxsize=32)
def _read_gtf_cached(path: str) -> pd.DataFrame:
    import pyranges1 as pr

    gr = pr.read_gtf(path)
    if hasattr(gr, "df"):
        return gr.df.copy()
    if hasattr(gr, "as_df"):
        return gr.as_df().copy()
    return pd.DataFrame(gr)


def _is_bigbed(path: str) -> bool:
    lowered = path.lower()
    if lowered.endswith(".bb") or lowered.endswith(".bigbed"):
        return True

    try:
        with pybigtools.open(path) as handle:
            is_bigbed = getattr(handle, "is_bigbed", None)
            return bool(is_bigbed() if callable(is_bigbed) else is_bigbed)
    except Exception:
        return False


def read_bed_regions(path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Read BED-like regions using indexed BigBed or cached pyranges1 BED parsing."""
    if _is_bigbed(path):
        with pybigtools.open(path) as bigbed:
            records = list(bigbed.records(chrom, start, end))
        return intervals_to_dataframe(records, chrom)

    bed_df = _read_bed_cached(path)
    return _slice_region(bed_df, chrom, start, end)


def read_gtf_regions(path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Read GTF regions using cached pyranges1 parsing."""
    gtf_df = _read_gtf_cached(path)
    return _slice_region(gtf_df, chrom, start, end)
