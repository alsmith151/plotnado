"""
Scaling utilities for genomic tracks.
"""

from typing import List, Literal
import numpy as np
import pandas as pd
from .base import GenomicRegion, Track


class Autoscaler:
    """
    Autoscale the data from multiple tracks to a single scale.
    """

    def __init__(
        self,
        tracks: List[Track],
        gr: GenomicRegion,
    ):
        self.tracks = tracks
        self.gr = gr

    @property
    def data(self) -> np.ndarray:
        """
        Get the data from all tracks for the specified region.
        """
        _data = []
        for t in self.tracks:
            data = t.fetch_data(self.gr)
            if isinstance(data, pd.DataFrame):
                # Assume the signal value is in the 'value' column
                if "value" in data.columns:
                    values = data["value"].values
                else:
                    # Fallback to last column
                    values = data.values[:, -1]
                _data.append(values)
            elif isinstance(data, np.ndarray):
                _data.append(data.flatten())

        if not _data:
            return np.array([0])

        return np.concatenate(_data, axis=0)

    @property
    def max_value(self) -> float:
        """Maximum value across all tracks."""
        d = self.data
        return float(np.nanmax(d)) if d.size > 0 else 1.0

    @property
    def min_value(self) -> float:
        """Minimum value across all tracks (bounded at 0 if all positive)."""
        d = self.data
        if d.size == 0:
            return 0.0
        min_val = np.nanmin(d)
        return float(min(0, min_val)) if min_val >= 0 else float(min_val)

    @property
    def mean_value(self) -> float:
        """Mean value across all tracks."""
        d = self.data
        return float(np.nanmean(d)) if d.size > 0 else 0.5


class Scaler:
    """
    Calculate scaling factors for tracks based on different methods.
    """

    def __init__(
        self,
        tracks: List[Track],
        gr: GenomicRegion,
        method: Literal["max", "mean", "total"] = "mean",
    ):
        self.tracks = tracks
        self.gr = gr
        self.method = method

    @property
    def data(self) -> List[np.ndarray]:
        """Fetch data arrays for each track."""
        _data = []
        for t in self.tracks:
            data = t.fetch_data(self.gr)
            if isinstance(data, pd.DataFrame):
                values = (
                    data["value"].values
                    if "value" in data.columns
                    else data.values[:, -1]
                )
                _data.append(values)
            elif isinstance(data, np.ndarray):
                _data.append(data.flatten())
        return _data

    @property
    def scaling_factors(self) -> np.ndarray:
        """Calculate scaling factors relative to the maximum across tracks."""
        data_list = self.data
        if not data_list:
            return np.array([])

        if self.method == "max":
            arr = [np.nanmax(d) if d.size > 0 else 0 for d in data_list]
        elif self.method == "mean":
            arr = [np.nanmean(d) if d.size > 0 else 0 for d in data_list]
        elif self.method == "total":
            arr = [np.nansum(d) if d.size > 0 else 0 for d in data_list]
        else:
            arr = [1.0] * len(data_list)

        max_val = np.max(arr)
        if max_val == 0:
            return np.ones(len(arr))

        return np.array(arr) / max_val
