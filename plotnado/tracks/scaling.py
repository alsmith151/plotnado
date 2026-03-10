"""
Scaling utilities for genomic tracks.
"""

import numpy as np
import pandas as pd
from .base import GenomicRegion, Track
from .enums import ScalingMethod


class Autoscaler:
    """
    Autoscale the data from multiple tracks to a single scale.
    """

    def __init__(
        self,
        tracks: list[Track],
        gr: GenomicRegion,
    ):
        self.tracks = tracks
        self.gr = gr

    @staticmethod
    def _supports_autoscale(track: Track) -> bool:
        return hasattr(track, "aesthetics") and (
            hasattr(track.aesthetics, "min_value") or hasattr(track.aesthetics, "max_value")
        )

    @staticmethod
    def _extract_numeric_values(data: pd.DataFrame | np.ndarray) -> np.ndarray:
        if isinstance(data, pd.DataFrame):
            if "value" in data.columns:
                series = pd.to_numeric(data["value"], errors="coerce")
                return series.dropna().to_numpy(dtype=float)
            return np.array([], dtype=float)

        flattened = np.asarray(data).ravel()
        if flattened.size == 0:
            return np.array([], dtype=float)
        numeric = pd.to_numeric(pd.Series(flattened), errors="coerce").dropna()
        return numeric.to_numpy(dtype=float)

    @property
    def data(self) -> np.ndarray:
        """
        Get the data from all tracks for the specified region.
        """
        _data = []
        for t in self.tracks:
            if not self._supports_autoscale(t):
                continue
            data = t.fetch_data(self.gr)
            if isinstance(data, (pd.DataFrame, np.ndarray)):
                values = self._extract_numeric_values(data)
                if values.size > 0:
                    _data.append(values)

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

    def apply(self) -> None:
        """Apply global min/max scale to compatible tracks."""
        min_value = self.min_value
        max_value = self.max_value
        if min_value == max_value:
            max_value = min_value + 1

        for track in self.tracks:
            if not self._supports_autoscale(track):
                continue
            if hasattr(track.aesthetics, "min_value"):
                track.aesthetics.min_value = min_value
            if hasattr(track.aesthetics, "max_value"):
                track.aesthetics.max_value = max_value


class Scaler:
    """
    Calculate scaling factors for tracks based on different methods.
    """

    def __init__(
        self,
        tracks: list[Track],
        gr: GenomicRegion,
        method: ScalingMethod = ScalingMethod.MEAN,
    ):
        self.tracks = tracks
        self.gr = gr
        self.method = method

    @property
    def data(self) -> list[np.ndarray]:
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

        if self.method == ScalingMethod.MAX:
            arr = [np.nanmax(d) if d.size > 0 else 0 for d in data_list]
        elif self.method == ScalingMethod.MEAN:
            arr = [np.nanmean(d) if d.size > 0 else 0 for d in data_list]
        elif self.method == ScalingMethod.TOTAL:
            arr = [np.nansum(d) if d.size > 0 else 0 for d in data_list]
        else:
            arr = [1.0] * len(data_list)

        max_val = np.max(arr)
        if max_val == 0:
            return np.ones(len(arr))

        return np.array(arr) / max_val
