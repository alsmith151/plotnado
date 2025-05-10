import functools
import math
import os
import pathlib
from enum import Enum
from typing import List, Literal

import coolbox.api as cb
import cooler.api as cooler
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import pyranges as pr
import tqdm
from coolbox.api import GenomeRange
from coolbox.core.track import Track
from coolbox.utilities import get_coverage_stack, get_feature_stack
from matplotlib import cm, transforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from pybedtools import BedTool
from typing import Optional, Union, List, Dict

from plotnado.api.utils import (
    get_human_readable_number_of_bp,
    interval_to_pyranges,
)

from .patches import plot_label, plot_text_range


class Autoscaler:
    """
    Autoscale the data from multiple tracks to a single scale.

    Args:
        tracks (List[cb.Track]): List of tracks to autoscale
        gr (cb.GenomeRange): Genome range to fetch the data
        gr2 (cb.GenomeRange, optional): Second genome range for 2D data. Defaults to None.

    """

    def __init__(
        self,
        tracks: List[cb.Track],
        gr: cb.GenomeRange,
        gr2: cb.GenomeRange = None,
    ):
        self.tracks = tracks
        self.gr = gr
        self.gr2 = gr2

        from .track_wrapper import MATRIX_TRACKS, BIGWIG_TRACKS, TrackWrapper

        assert len(tracks) > 0, "No tracks to autoscale"
        assert all(
            isinstance(t, (cb.Track, TrackWrapper)) for t in tracks
        ), "All tracks must be of type cb.Track"
        assert all(
            type(t) in MATRIX_TRACKS + BIGWIG_TRACKS for t in tracks
        ), "All tracks must be of tracks that produce numerical data (MatrixCapcruncher, BigWig, etc)"

    @property
    def data(self):
        """
        Get the data from all tracks for the specified region
        """
        _data = [t.fetch_data(gr=self.gr, gr2=self.gr2) for t in self.tracks]

        if isinstance(
            _data[0], pd.DataFrame
        ):  # If the data is a DataFrame, we need to extract the values from the last column
            _data = [d.values[:, -1] for d in _data]
        elif isinstance(_data[0], np.ndarray):
            pass

        return np.concatenate(_data, axis=0)

    @property
    def max_value(self):
        return np.nanmax(self.data, axis=0)

    @property
    def min_value(self):
        return min(0, np.nanmin(self.data, axis=0))

    @property
    def mean_value(self):
        return np.nanmean(self.data, axis=0)


class Scaler:
    def __init__(
        self,
        tracks: List[cb.Track],
        gr: cb.GenomeRange,
        gr2: cb.GenomeRange = None,
        method: Literal["max", "mean", "total"] = "mean",
    ):
        self.tracks = tracks
        self.gr = gr
        self.gr2 = gr2
        self.method = method

        from .track_wrapper import MATRIX_TRACKS, BIGWIG_TRACKS, TrackWrapper

        assert len(tracks) > 0, "No tracks to autoscale"
        assert all(
            isinstance(t, (cb.Track, TrackWrapper)) for t in tracks
        ), "All tracks must be of type cb.Track"
        assert all(
            type(t) in MATRIX_TRACKS + BIGWIG_TRACKS for t in tracks
        ), "All tracks must be of tracks that produce numerical data (MatrixCapcruncher, BigWig, etc)"

    @property
    def data(self) -> List[np.ndarray]:
        """
        Get the data from all tracks for the specified region
        """
        _data = [t.fetch_data(gr=self.gr, gr2=self.gr2) for t in self.tracks]

        if isinstance(
            _data[0], pd.DataFrame
        ):  # If the data is a DataFrame, we need to extract the values from the last column
            _data = [d.values[:, -1] for d in _data]
        elif isinstance(_data[0], np.ndarray):
            pass

        return _data

    @property
    def scaling_factors(self) -> np.ndarray:
        if self.method == "max":
            arr = [np.nanmax(d, axis=0) for d in self.data]
        elif self.method == "mean":
            arr = [np.nanmean(d, axis=0) for d in self.data]
        elif self.method == "total":
            arr = [np.nansum(d, axis=0) for d in self.data]

        return np.array(arr) / np.max(arr)


class MatrixCapcruncher(cb.Cool):
    """
    Matrix track designed to plot CapCruncher derived matrices.

    Args:
        file (os.PathLike): Path to the cooler file
        binsize (int, optional): Binsize of the matrix. Defaults to 5000.
        viewpoint (str): Viewpoint to plot the matrix from
        remove_viewpoint (bool, optional): Remove the viewpoint from the matrix. Defaults to False.
        **kwargs: Additional arguments to pass to the track
    """

    def __init__(
        self,
        file: os.PathLike,
        binsize: 5000,
        viewpoint: str,
        remove_viewpoint=False,
        **kwargs,
    ):
        self.binsize = binsize
        self.viewpoint = viewpoint
        self.remove_viewpoint = remove_viewpoint
        self.properties = dict()
        self.properties.update(kwargs)
        self.properties["name"] = f"CCMatrix.{self.properties.get('title')}"
        super(MatrixCapcruncher, self).__init__(file, **kwargs)
        # Need to override the coolbox default if we need a cmap to be set
        self.properties["color"] = kwargs.get("color", self.properties["color"])

        # Override the defaults
        self.properties["balance"] = "no"

        if not self._cooler_store_has_binsize:
            raise ValueError(
                f"Viewpoint {viewpoint} or resolution {binsize} not found in supplied file."
            )

        self.cooler = cooler.Cooler(f"{file}::{viewpoint}/resolutions/{binsize}")
        self.capture_bins = self.cooler.info["metadata"]["viewpoint_bins"]

    def _cooler_store_has_binsize(self):
        clrs = cooler.fileops.list_coolers(self.file)
        expected_path = f"{self.viewpoint}/resolutions/{self.binsize}"

        if expected_path in clrs:
            return True

    def get_matrix(self, coordinates, field="count"):
        matrix = self.cooler.matrix(field=field, balance=False).fetch(coordinates)

        offset = self.cooler.offset(coordinates)
        capture_bins = [(bin - offset) for bin in self.capture_bins]

        if self.remove_viewpoint:
            matrix[capture_bins, :] = 0
            matrix[:, capture_bins] = 0

        return matrix

    def get_matrix_normalised(
        self, coordinates, normalization_method=None, **normalisation_kwargs
    ):
        methods_stored = {
            "n_interactions": "count_n_interactions_norm",
            "n_rf_n_interactions": "count_n_rf_n_interactions_norm",
        }

        if normalization_method == "raw":
            matrix_normalised = self.get_matrix(coordinates)

        elif normalization_method in methods_stored:
            matrix_normalised = self.get_matrix(
                coordinates, field=methods_stored[normalization_method]
            )

        elif normalization_method == "ice":
            import iced

            matrix = self.get_matrix(coordinates)
            matrix = np.nan_to_num(matrix)
            # matrix = iced.filter.filter_low_counts(matrix, percentage=0.04)
            matrix_normalised = iced.normalization.ICE_normalization(
                matrix, **normalisation_kwargs
            )  # Get iced matrix

        elif normalization_method == "icen_cis":
            import iced

            matrix = self.get_matrix(coordinates)
            matrix = np.nan_to_num(matrix)
            matrix_ice = iced.normalization.ICE_normalization(
                matrix, **normalisation_kwargs
            )  # Get iced matrix
            matrix_normalised = (
                matrix_ice
                / int(self.cooler.info["metadata"]["n_cis_interactions"])
                * 1e6
            )  # Correct for number of interactions * 1e6

        elif normalization_method == "icen_scale":
            matrix = self.get_matrix(coordinates)
            matrix = np.nan_to_num(matrix)
            matrix_ice = iced.normalization.ICE_normalization(
                matrix, **normalisation_kwargs
            )  # Get iced matrix
            matrix_normalised = matrix_ice / self.properties["scaling_factor"]

        else:
            raise ValueError(
                f'Incorrect normalisation specified choose from: {" ".join(["raw", *methods_stored.keys(),"ice", "icen_cis", "icen_scale"])}'
            )

        return matrix_normalised

    def fetch_data(
        self, gr: cb.GenomeRange, gr2: cb.GenomeRange = None, **kwargs
    ) -> np.ndarray:
        norm = self.properties.get("normalization", "raw")
        matrix = self.get_matrix_normalised(
            f"{gr.chrom}:{gr.start}-{gr.end}", normalization_method=norm, **kwargs
        )
        return self.fill_zero_nan(matrix)

    def plot_matrix(self, gr: GenomeRange, gr2: GenomeRange = None):
        # Code taken and adapted from coolbox
        gr = GenomeRange(gr)

        if "JuiceBox" in self.properties["color"]:
            cmap = MatrixCapcruncher.get_juicebox_cmaps()[self.properties["color"]]
        else:
            cmap = cm.get_cmap(self.properties["color"])

        lowest = cmap(0)
        cmap.set_bad(lowest)
        cmap.set_under(lowest)

        ax = self.ax
        arr = self.matrix
        c_min, c_max = self.matrix_val_range

        if self.properties["max_value"] == "auto":
            matrix_triu = np.triu(self.matrix)
            c_max = np.percentile(matrix_triu, 98)

        if gr2 is None and self.style == self.STYLE_TRIANGULAR:
            # triangular style
            scale_r = 1 / math.sqrt(2)
            r_len = gr.end - gr.start
            # Rotate image using Affine2D, reference:
            #     https://stackoverflow.com/a/50920567/8500469

            tr = (
                transforms.Affine2D()
                .translate(-gr.start, -gr.start)
                .rotate_deg_around(0, 0, 45)
                .scale(scale_r)
                .translate(gr.start + r_len / 2, -r_len / 2)
            )

            img = ax.matshow(
                arr,
                cmap=cmap,
                transform=tr + ax.transData,
                extent=(gr.start, gr.end, gr.start, gr.end),
                aspect="auto",
                interpolation="none",
            )

        elif gr2 is None and self.style == self.STYLE_WINDOW:
            # window style
            # exist in HicMatBase
            fgr = self.fetched_gr
            scale_factor = fgr.length / gr.length
            scale_r = scale_factor / math.sqrt(2)
            length_dialog = gr.length * scale_factor
            delta_x = length_dialog * (gr.start - fgr.start) / fgr.length
            delta_x = length_dialog / 2 - delta_x
            tr = (
                transforms.Affine2D()
                .translate(-gr.start, -gr.start)
                .rotate_deg_around(0, 0, 45)
                .scale(scale_r)
                .translate(gr.start + delta_x, -fgr.length / 2)
            )
            img = ax.matshow(
                arr,
                cmap=cmap,
                transform=tr + ax.transData,
                extent=(gr.start, gr.end, gr.start, gr.end),
                aspect="auto",
            )
        else:
            if gr2 is None:
                gr2 = gr
            # matrix style
            img = ax.matshow(
                arr,
                cmap=cmap,
                extent=(gr.start, gr.end, gr2.end, gr2.start),
                aspect="auto",
            )

        if self.norm == "log":
            img.set_norm(colors.LogNorm(vmin=c_min, vmax=c_max))
        else:
            img.set_norm(colors.Normalize(vmin=c_min, vmax=c_max))

        return img

    @staticmethod
    def get_juicebox_cmaps():
        JuiceBoxLikeColor = LinearSegmentedColormap.from_list(
            "interaction", ["#FFFFFF", "#FFDFDF", "#FF7575", "#FF2626", "#F70000"]
        )
        JuiceBoxLikeColor.set_bad("white")
        JuiceBoxLikeColor.set_under("white")
        JuiceBoxLikeColor2 = LinearSegmentedColormap.from_list(
            "interaction", ["#FFFFFF", "#FFDFAF", "#FF7555", "#FF2600", "#F70000"]
        )
        JuiceBoxLikeColor2.set_bad("white")
        JuiceBoxLikeColor2.set_under("white")

        return {
            "JuiceBoxLike": JuiceBoxLikeColor,
            "JuiceBoxLike2": JuiceBoxLikeColor2,
        }


class BigwigFragment(cb.BigWig):
    """
    Subclass of BigWig that plots fragments instead of lines. Provides a more accurate representation of the data.

    Args:
        file (os.PathLike): Path to the bigwig file
        **kwargs: Additional arguments to pass to the track
    """

    def __init__(self, file, **kwargs):
        self.file = file
        self.coverages = []

        super(BigwigFragment, self).__init__(file, **kwargs)

    def fetch_data(self, gr, **kwargs):
        if not self.properties["style"] == "fragment":
            data = super(BigwigFragment, self).fetch_data(gr, **kwargs)
        else:
            data = self.bw.fetch_intervals(gr.chrom, gr.start, gr.end)

        return data

    def plot_fragments(self, ax, gr, **kwargs):
        data = self.fetch_data(gr, **kwargs)
        _alpha = self.properties.get("alpha", 1.0)
        _threshold = self.properties.get("threshold", 0)
        _offset = gr.start
        bp_proportion = 1 / (data["end"].max() - data["start"].min())

        for row in data.itertuples():
            pg = Polygon(
                [
                    (row.start, 0),
                    (row.start, row.value),
                    (row.end, row.value),
                    (row.end, 0),
                ],
                color=self.properties["color"],
            )
            ax.add_patch(pg)

        ax.set_ylim(0, data["value"].max())
        ymin, ymax = self.adjust_plot(ax, gr)

        if self.properties.get("show_data_range") == "yes" and not self.properties.get(
            "is_subtrack"
        ):
            self.plot_data_range(
                ax, ymin, ymax, self.properties["data_range_style"], gr
            )

        if self.properties.get("label_on_track") and not self.properties.get(
            "is_subtrack"
        ):
            self.plot_label()

    def plot(self, ax, gr, **kwargs):
        if not self.properties["style"] == "fragment":
            super(BigwigFragment, self).plot(ax, gr, **kwargs)
        else:
            self.plot_fragments(ax, gr, **kwargs)


class BigwigFragmentCollection(Track):
    """
    Aggregates multiple bigwig files into a single track. Plots mean and standard error of the mean.


    Args:
        file (list): List of paths to the bigwig files
        exclusions (str, optional): Path to a bed file containing regions to exclude. Defaults to None.
        **kwargs: Additional arguments to pass to the track

    """

    DEFAULT_PROPERTIES = {
        "style": "line",
        "fmt": "-",
        "line_width": 2.0,
        "size": 10,
        "color": "#a6cee3",
        "threshold_color": "#ff9c9c",
        "threshold": "inf",
        "cmap": "bwr",
        "orientation": None,
        "data_range_style": "y-axis",
        "min_value": "auto",
        "max_value": "auto",
    }

    def __init__(self, file: list, exclusions: str = None, **kwargs):
        self.file_names = file
        self.exclusions = exclusions
        self.bws = (
            [cb.BigWig(str(fn)) for fn in file]
            if not isinstance(file[0], cb.BigWig)
            else file
        )
        self.properties = {"files": self.file_names}
        self.properties.update(BigwigFragmentCollection.DEFAULT_PROPERTIES.copy())
        self.properties.update(kwargs)
        self.properties["name"] = f"BigWigCollection.{self.properties.get('title')}"
        super(BigwigFragmentCollection, self).__init__(**self.properties)

        self.coverages = []

        # load features from global feature stack
        features_stack = get_feature_stack()
        for features in features_stack:
            self.properties.update(features.properties)

        # load coverages from global coverages stack
        coverage_stack = get_coverage_stack()
        for coverage in coverage_stack:
            self.coverages.append(coverage)

    @property
    def subtracks(self):
        return self.bws

    def fetch_data(self, gr, **kwargs):
        datasets = [
            bw.bw.fetch_intervals(gr.chrom, gr.start, gr.end)
            .set_index(["chrom", "start", "end"])
            .rename(columns={"value": os.path.basename(bw.properties["file"])})
            for bw in self.bws
        ]
        df = datasets[0].join(datasets[1:])
        df_summary = df.assign(mean=df.mean(axis=1), sem=df.sem(axis=1)).reset_index()

        intervals_to_bp = []
        for interval in df_summary.itertuples():
            interval_len = interval.end - interval.start

            interval_positions = np.arange(interval_len) + interval.start
            scores_mean = np.repeat(interval.mean, interval_len)
            scores_sem = np.repeat(interval.sem, interval_len)

            intervals_to_bp.append(
                np.vstack([interval_positions, scores_mean, scores_sem]).T
            )

        df_intervals = pd.DataFrame(
            np.concatenate(intervals_to_bp), columns=["bp", "mean", "sem"]
        )

        if self.exclusions:
            df_intervals = pd.concat(
                [df_intervals, self.fetch_exluded_regions(gr)]
            ).sort_values("bp")

        if self.properties.get("smooth_window") or self.properties.get("smooth"):
            from scipy.signal import savgol_filter

            df_intervals["mean_smoothed"] = savgol_filter(
                df_intervals["mean"],
                window_length=self.properties.get("smooth_window", 1001),
                polyorder=self.properties.get("polyorder", 1),
            )

        return df_intervals

    def fetch_exluded_regions(self, gr):
        excluded_tabix = BedTool(self.exclusions).tabix(force=True)
        df_excluded = excluded_tabix.tabix_intervals(
            f"{gr.chrom}:{gr.start}-{gr.end}"
        ).to_dataframe()

        intervals_to_bp = []
        for interval in df_excluded.itertuples():
            interval_len = interval.end - interval.start

            interval_positions = np.arange(interval_len) + interval.start
            scores_nan = np.repeat(np.nan, interval_len)
            intervals_to_bp.append(interval_positions)

        df_intervals = pd.Series(np.concatenate(intervals_to_bp)).to_frame("bp")
        df_intervals["mean"] = np.nan
        df_intervals["sem"] = np.nan

        return df_intervals

    def plot(self, ax, gr, **kwargs):
        data = self.fetch_data(gr, **kwargs)

        line_width = self.properties.get("line_width", 1)
        color = self.properties.get("color", "blue")
        alpha = self.properties.get("alpha", 0.2)
        downsample = self.properties.get(
            "downsample", 10
        )  # downsample the data to make it faster to plot

        if downsample:
            rows = np.arange(0, data.shape[0], downsample)
            data = data.iloc[rows]

        if self.properties.get("smooth_window"):
            scores = data["mean_smoothed"]
        else:
            scores = data["mean"]

        ax.fill_between(
            data["bp"],
            scores - data["sem"],
            scores + data["sem"],
            alpha=alpha,
            color=color,
            zorder=0,
        )

        ax.plot(
            data["bp"],
            scores,
            color=color,
            zorder=1,
            linewidth=line_width,
        )

        min_val = self.properties.get("min_value")
        max_val = self.properties.get("max_value")

        ymin = round(scores.min()) if min_val == "auto" else min_val
        ymin = min(0, ymin)

        ymax = round(scores.max() + data["sem"].max()) if max_val == "auto" else max_val

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(ymin, ymax)

        if not self.properties.get("is_subtrack"):
            self.plot_data_range(
                ax, ymin, ymax, self.properties["data_range_style"], gr
            )

        if self.properties.get("label_on_track") and not self.properties.get(
            "is_subtrack"
        ):
            self.plot_label()

    def plot_data_range(self, ax, ymin, ymax, data_range_style, gr: cb.GenomeRange):
        if data_range_style == "text":
            self._plot_text_range(ax, ymin, ymax, gr)
        else:  # 'y-axis' style
            try:
                y_ax = self.y_ax
                self.plot_yaxis_range(ax, y_ax)
            except AttributeError:
                self.plot_data_range(ax, ymin, ymax, "text", gr)

    def _plot_text_range(self, ax, ymin, ymax, gr: cb.GenomeRange):
        plot_text_range(self, ax, ymin, ymax, gr)

    def plot_yaxis_range(self, plot_axis, y_ax):
        # """
        # Plot the scale of the y axis with respect to the plot_axis
        # plot something that looks like this:
        # ymax ┐
        #      │
        #      │
        # ymin ┘
        # Parameters
        # ----------
        # plot_axis : matplotlib.axes.Axes
        #     Main plot axis.
        # y_ax : matplotlib.axes.Axes
        #     Axis to use to plot the scale
        # """

        if (
            "show_data_range" in self.properties
            and self.properties["show_data_range"] == "no"
        ):
            return

        def value_to_str(value):
            if value % 1 == 0:
                return str(int(value))
            else:
                return f"{value:.4f}" if value < 0.01 else f"{value:.2f}"

        ymin, ymax = plot_axis.get_ylim()

        ymax_str = value_to_str(ymax)
        ymin_str = value_to_str(ymin)
        x_pos = [0, 0.5, 0.5, 0]
        y_pos = [0.01, 0.01, 0.99, 0.99]
        y_ax.plot(x_pos, y_pos, color="black", linewidth=1, transform=y_ax.transAxes)
        y_ax.text(
            -0.2,
            -0.01,
            ymin_str,
            verticalalignment="bottom",
            horizontalalignment="right",
            transform=y_ax.transAxes,
        )
        y_ax.text(
            -0.2,
            1,
            ymax_str,
            verticalalignment="top",
            horizontalalignment="right",
            transform=y_ax.transAxes,
        )
        y_ax.patch.set_visible(False)

    def plot_text_range(self, ax, ymin, ymax, gr: cb.GenomeRange):
        ydelta = ymax - ymin

        # set min max
        def format_lim(lim):
            return int(lim) if float(lim) % 1 == 0 else f"{lim:.2f}"

        ymax_print = format_lim(ymax)
        ymin_print = format_lim(ymin)
        small_x = 0.01 * gr.length
        # by default show the data range
        ax.text(
            gr.start - small_x,
            ymax - ydelta * 0.2,
            f"[ {ymin_print} ~ {ymax_print} ]",
            horizontalalignment="left",
            verticalalignment="top",
        )


class BigwigFragmentCollectionOverlay(Track):
    """
    Overlay multiple bigwig collections on top of each other.

    Args:
        collections (List[FragmentBigwigCollection]): List of FragmentBigwigCollection objects
        **kwargs: Additional arguments to pass to the track
    """

    def __init__(self, collections: List[BigwigFragmentCollection], **kwargs):
        # Initialization
        self.collections = collections

        for track in self.collections:
            track.properties["show_data_range"] = "no"
            track.properties["data_range_style"] = "no"
            track.properties["label_on_track"] = "no"
            track.properties["is_subtrack"] = True

        self.properties = {"collections": collections}
        self.properties.update(kwargs)

    @property
    def subtracks(self):
        return self.collections

    def plot(self, ax, genome_range, **kwargs):
        # Method for plot the track
        # within the genome range "chrom:start-end". `genome_range` is a `GenomeRange` object with properties: chrom, start, end
        # Draw in the pass-in `ax` Axes
        for collection in self.collections:
            collection.plot(ax, genome_range, **kwargs)

        self.plot_label()

    def fetch_data(self, genome_range, **kwargs):
        # Method for fetch the data within genome_range
        # genome_range is a `coolbox.utilities.GenomeRange` object
        for collection in self.collections:
            collection.fetch_data(genome_range, **kwargs)


def plot_bigwig_scaled(track, ax, gr, scale_factor, **kwargs):
    genome_range = gr

    data = track.fetch_plot_data(gr, **kwargs)
    if isinstance(data, pd.DataFrame):
        # ['pos', 'score']
        if "pos" in data:
            data = data["pos"].to_numpy(), data["score"].to_numpy()
        # ['score']
        else:
            data = np.linspace(gr.start, gr.end, len(data)), data["score"].to_numpy()
    elif isinstance(data, np.ndarray):
        # 1d or 2d ndarray
        if len(data.shape) == 1:
            data = np.linspace(gr.start, gr.end, len(data)), data
        elif len(data.shape) == 2:
            data = np.linspace(gr.start, gr.end, data.shape[1]), data
        else:
            raise ValueError("The ndarray must be in 1d or 2d format")
    elif not isinstance(data, tuple):
        # not (indexes, values)
        raise ValueError(f"Data format not supported.")

    indexes, values = data
    values = values / scale_factor
    track.plot_hist(ax, gr, indexes, values)
    track.plot_label()


class BigwigOverlay(Track):
    """
    Overlay multiple bigwig files on top of each other.

    Args:
        files (List[os.PathLike]): List of paths to the bigwig files
        **kwargs: Additional arguments to pass to the track
    """

    def __init__(self, collection: List[Union[cb.BigWig, str]], **kwargs):
        _collection = []
        for obj in collection:
            if isinstance(obj, str):
                assert os.path.exists(obj) or kwargs.get(
                    "ignore_file_validation"
                ), f"File {obj} does not exist"
                _collection.append(cb.BigWig(obj))
            elif isinstance(obj, cb.Track):
                _collection.append(obj)
            else:
                raise ValueError(
                    f"Object {obj} is not a valid BigWig or path to a BigWig file"
                )

        self.coverages = []
        self.collection = _collection
        self.properties = dict()
        self.properties.update(
            {
                "min_value": "auto",
                "max_value": "auto",
                "scale": False,
                "scale_n_cis": False,
            }
        )
        self.properties.update(kwargs)
        self.properties["name"] = f"BigWigOverlay.{self.properties.get('title')}"

    @property
    def subtracks(self):
        return self.collection

    def fetch_data(
        self, gr: GenomeRange, scale: bool = False, scale_n_cis: bool = False, **kwargs
    ):
        data = []

        if not scale:
            for bw in self.collection:
                data.append(bw.fetch_data(gr, **kwargs))

        elif scale and not scale_n_cis:
            scaler = Scaler(tracks=self.collection, gr=gr)
            scaling_factors = scaler.scaling_factors

            for bw, scaling_factor in zip(self.collection, scaling_factors):
                data.append(bw.fetch_data(gr, **kwargs) / scaling_factor)

        elif scale and scale_n_cis:
            import pyBigWig

            chrom = gr.chrom
            if isinstance(self.collection[0], cb.BigWig):
                bw_file = self.collection[0].properties["file"]
            else:
                bw_file = self.collection[0]

            with pyBigWig.open(bw_file) as bw:
                chrom_length = bw.chroms()[chrom]

            gr_chrom = GenomeRange(chrom, 0, chrom_length)
            scaler = Scaler(tracks=self.collection, gr=gr_chrom, method="total")
            scaling_factors = scaler.scaling_factors

            for bw, scaling_factor in zip(self.collection, scaling_factors):
                data.append(bw.fetch_data(gr, **kwargs) / scaling_factor)

        return data

    def _get_scaling_factors(
        self, gr: GenomeRange, method: Literal["n_cis", "mean", "max", "min"]
    ):
        if method == "n_cis":
            import pyBigWig

            chrom = gr.chrom
            if isinstance(self.collection[0], cb.BigWig):
                bw_file = self.collection[0].properties["file"]
            else:
                bw_file = self.collection[0]

            with pyBigWig.open(bw_file) as bw:
                chrom_length = bw.chroms()[chrom]

            gr_chrom = GenomeRange(chrom, 0, chrom_length)
            scaler = Scaler(tracks=self.collection, gr=gr_chrom, method="total")
            scaling_factors = scaler.scaling_factors
        else:
            scaler = Scaler(tracks=self.collection, gr=gr, method=method)
            scaling_factors = scaler.scaling_factors

        return scaling_factors

    def plot(self, ax, gr: GenomeRange, **kwargs):
        scaler = Autoscaler(tracks=self.collection, gr=gr)
        min_value = (
            scaler.min_value
            if self.properties.get("min_value") == "auto"
            else self.properties.get("min_value", 0)
        )
        max_value = (
            scaler.max_value
            if self.properties.get("max_value") == "auto"
            else self.properties.get("max_value")
        )

        if self.properties.get("scale_method"):
            scaling_factors = self._get_scaling_factors(
                gr, self.properties.get("scale_method")
            )
            for ii, (bw, scaling_factor) in enumerate(
                zip(self.collection, scaling_factors)
            ):
                bw.properties["min_value"] = min_value
                bw.properties["max_value"] = max_value
                bw.properties["show_data_range"] = "no"
                bw.properties["data_range_style"] = "no"
                bw.properties["scaling_factor"] = scaling_factor
                plot_bigwig_scaled(bw, ax, gr, scaling_factor, **kwargs)

        else:
            for ii, bw in enumerate(self.collection):
                bw.properties["show_data_range"] = "no"
                bw.properties["data_range_style"] = "no"
                bw.properties["min_value"] = min_value
                bw.properties["max_value"] = max_value
                bw.properties["is_subtrack"] = True
                bw.plot(
                    ax,
                    gr,
                    show_data_range=False,
                    **kwargs,
                )

        self._plot_text_range(ax, min_value, max_value, gr)
        self.plot_label()

    def _plot_text_range(self, ax, ymin, ymax, gr: cb.GenomeRange):
        plot_text_range(self, ax, ymin, ymax, gr)

    def plot_label(self):
        if self.properties.get("label_on_track") not in [
            "True",
            "yes",
            "T",
            "Y",
            "1",
            True,
            1,
        ]:
            if (
                hasattr(self, "label_ax")
                and self.label_ax is not None
                and "title" in self.properties
            ):
                self.label_ax.text(
                    0.15,
                    0.5,
                    self.properties["title"],
                    horizontalalignment="left",
                    size="large",
                    verticalalignment="center",
                )
            elif hasattr(self, "label_ax") and self.properties.get("label_subtracks"):
                # Plot a patch with the color of the subtracks and the name of the subtrack
                for i, bw in enumerate(self.collection):
                    self.label_ax.add_patch(
                        Polygon(
                            [
                                (0.1, 0.5 - (i * 0.1)),
                                (0.1, 0.6 - (i * 0.1)),
                                (0.2, 0.6 - (i * 0.1)),
                                (0.2, 0.5 - (i * 0.1)),
                            ],
                            color=bw.properties.get("color", "black"),
                        )
                    )

                    self.label_ax.text(
                        0.25,
                        0.55 - (i * 0.1),
                        bw.properties.get("title", f"Subtrack {i}"),
                        horizontalalignment="left",
                        verticalalignment="center",
                    )

    def adjust_plot(self, ax, gr: GenomeRange):
        ax.set_xlim(gr.start, gr.end)
        ymin, ymax = ax.get_ylim()
        if "max_value" in self.properties and self.properties["max_value"] != "auto":
            ymax = self.properties["max_value"]
        if "min_value" in self.properties and self.properties["min_value"] != "auto":
            ymin = self.properties["min_value"]

        if (
            "orientation" in self.properties
            and self.properties["orientation"] == "inverted"
        ):
            ax.set_ylim(ymax, ymin)
        else:
            ax.set_ylim(ymin, ymax)
        return ymin, ymax


class BigwigSubtraction(cb.BigWig):
    """
    Compare two bigwig files by plotting their difference.

    Args:
        file1 (os.PathLike): Path to the first bigwig file
        file2 (os.PathLike): Path to the second bigwig file
        **kwargs: Additional arguments to pass to the track
    """

    def __init__(self, files: List[str], **kwargs):
        self.properties = dict()
        self.properties.update(kwargs)
        self.properties["name"] = f"BigWigSubtraction.{self.properties.get('title')}"
        super(BigwigSubtraction, self).__init__(files[0], **self.properties)

        self.bw1 = cb.BigWig(files[0])
        self.bw2 = cb.BigWig(files[1])

    @property
    def subtracks(self):
        return [self.bw1, self.bw2]

    def fetch_data(self, gr: GenomeRange, **kwargs):
        data1 = self.bw1.fetch_data(gr, **kwargs)
        data2 = self.bw2.fetch_data(gr, **kwargs)

        n_bins = len(range(gr.start, gr.end, self.properties.get("binsize", 10)))
        new_starts = np.linspace(gr.start, gr.end, n_bins)
        interpol_1 = np.interp(new_starts, data1.start, data1.score.values)
        interpol_2 = np.interp(new_starts, data2.start, data2.score.values)
        diff = interpol_1 - interpol_2

        return pd.DataFrame(
            {
                "chrom": gr.chrom,
                "start": new_starts,
                "end": pd.Series(new_starts).shift(-1),
                "score": diff,
            }
        ).dropna()

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        return self.fetch_data(gr, **kwargs)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        if not self.properties["style"] == "line":
            super(BigwigSubtraction, self).plot(ax, gr, **kwargs)
        else:
            data = self.fetch_data(gr, **kwargs)

            # Split the data into positive and negative values
            fmt = self.properties.get("fmt")
            line_width = self.properties.get("line_width", 1)
            color = self.properties.get("color", "black")
            alpha = self.properties.get("alpha", 1.0)
            threshold = float(self.properties.get("threshold", 0))
            threshold_color = self.properties.get("threshold_color")

            num_bins = self.properties.get("num_bins", 700)
            x = np.linspace(gr.start, gr.end, num_bins)
            y = np.interp(x, data["start"], data["score"])

            y_pos = y[y > threshold]
            y_neg = y[y < threshold]
            x_pos = x[y > threshold]
            x_neg = x[y < threshold]

            # PLot the positive values as a line
            ax.plot(
                x_pos,
                y_pos,
                color=color,
                zorder=1,
                alpha=alpha,
                linewidth=line_width,
                linestyle=fmt,
            )

            # Plot the negative values as a line
            ax.plot(
                x_neg,
                y_neg,
                color=threshold_color,
                zorder=1,
                alpha=alpha,
                linewidth=line_width,
                linestyle=fmt,
            )
            ymin, ymax = self.adjust_plot(ax, gr)
            # disable plot data range in coverage mode
            if (
                self.properties["data_range_style"] != "no"
                and "track" not in self.__dict__
            ):
                self.plot_data_range(
                    ax, ymin, ymax, self.properties["data_range_style"], gr
                )

            self.plot_label()

    def plot_label(self):
        if hasattr(self, "label_ax") and self.label_ax is not None:
            self.label_ax.text(
                0.15,
                0.5,
                self.properties["title"],
                horizontalalignment="left",
                size="large",
                verticalalignment="center",
            )


class BedSimple(cb.BED):
    """
    Simplified version of the Bed class that only plots the intervals as rectangles.

    Args:
        file (os.PathLike): Path to the bed file
    """

    def __init__(self, file: str, **kwargs):
        self.file = file
        self.properties = dict()
        self.properties["name"] = f"SimpleBed_{self.properties.get('title')}"
        self.properties.update(kwargs)

    def fetch_data(self, gr):
        bt = BedTool(self.file)

        try:
            bt_tabix = bt.tabix(force=True)
            intervals = bt_tabix.tabix_intervals(
                f"{gr.chrom}:{gr.start}-{gr.end}"
            )

        except OSError: # Handle the case where the bed file is not tabix indexed or the user does not have permission to write to the directory
            import tempfile
            with tempfile.NamedTemporaryFile() as tmp:
                bt.saveas(tmp.name)
                bt_tabix = BedTool(tmp.name).tabix(force=True)
                intervals = bt_tabix.tabix_intervals(
                    f"{gr.chrom}:{gr.start}-{gr.end}"
                )
        
        return intervals

    def plot(self, ax, gr, **kwargs):
        data = self.fetch_data(gr)
        y_midpoint = 0.5

        for interval in data.intervals:
            pg = Polygon(
                [
                    (interval.start, y_midpoint - 0.1),
                    (interval.start, y_midpoint + 0.1),
                    (interval.end, y_midpoint + 0.1),
                    (interval.end, y_midpoint - 0.1),
                ],
                color=self.properties.get("color", "black"),
            )

            if hasattr(interval, "name") and not self.properties.get("no_annotation"):
                interval_midpoint = interval.start + (
                    (interval.end - interval.start) / 2
                )
                ax.text(
                    interval_midpoint,
                    y_midpoint - 0.1,
                    interval.name,
                    ha="center",
                    va="center",
                )

            ax.add_patch(pg)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)

        if self.properties.get('title'):
            self.plot_label()

    def adjust_plot(self, ax, gr: GenomeRange):
        ax.set_xlim(gr.start, gr.end)
        ymin, ymax = ax.get_ylim()
        return ymin, ymax


class BedMemory(Track):
    """
    A bed track that is stored in memory as  pyranges object

    Args:
        pyranges (pr.PyRanges): A pyranges object
        **kwargs: Additional arguments to pass to the track

    """

    def __init__(self, pyranges: pr.PyRanges, **kwargs):
        self.file = None
        self.pyranges = pyranges
        self.properties = dict()
        self.properties["name"] = "MemmoryBed"
        self.properties.update(kwargs)

    def fetch_data(self, gr: GenomeRange):
        interval = interval_to_pyranges(gr)

        return self.pyranges.intersect(interval)

    def plot(self, ax, gr, **kwargs):
        data = self.fetch_data(gr)
        y_midpoint = 0.5

        for interval in data.as_df().itertuples():
            if hasattr(interval, "color"):
                color = interval.color
            else:
                color = self.properties.get("color", "black")

            pg = Polygon(
                [
                    (interval.Start, y_midpoint - 0.1),
                    (interval.Start, y_midpoint + 0.1),
                    (interval.End, y_midpoint + 0.1),
                    (interval.End, y_midpoint - 0.1),
                ],
                color=color,
            )

            if hasattr(interval, "Name") and not self.properties.get("no_annotation"):
                interval_midpoint = interval.Start + (
                    (interval.End - interval.Start) / 2
                )
                ax.text(
                    interval_midpoint,
                    y_midpoint - 0.1,
                    interval.Name,
                    ha="center",
                    va="center",
                )

            ax.add_patch(pg)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)

    def plot_label(self):
        if hasattr(self, "label_ax") and self.label_ax is not None:
            self.label_ax.text(
                0.15,
                0.5,
                self.properties["title"],
                horizontalalignment="left",
                size="large",
                verticalalignment="center",
            )


class ScaleBar(Track):
    """
    A scale bar that shows the length of the genomic region.
    """

    def __init__(self, **kwargs):
        self.properties = dict()
        self.properties["name"] = "Scale"
        self.properties.update(kwargs)
        super(ScaleBar, self).__init__()

    def fetch_data(self, **kwargs):
        pass

    def get_appropriate_scale(self, length):
        if length <= 1e3:
            scale = 1e2
        elif 1e3 < length < 1e4:
            scale = 1e3
        elif 1e4 < length < 1e5:
            scale = 1e4
        elif 1e5 < length < 1e6:
            scale = 1e5
        elif 1e6 < length < 1e7:
            scale = 1e6
        elif 1e7 < length < 1e8:
            scale = 1e8

        return scale

    def plot(self, ax, gr, **kwargs):
        position = self.properties.get("position", "left")
        y_midpoint = 0.5

        if self.properties.get("scale_distance"):
            scale_distance = self.properties["scale_distance"]

        else:
            scale_distance = self.get_appropriate_scale(gr.end - gr.start)

        # Determine x start and end based on position
        if position == "left":
            x0 = gr.start
            x1 = x0 + scale_distance
        elif position == "right":
            x0 = gr.end - scale_distance
            x1 = gr.end
        else:
            raise ValueError('Position can only be "left" or "right"')

        # Plot scale bar
        ax.plot([x0, x1], [y_midpoint, y_midpoint], color="black")
        ax.plot([x0, x0], [y_midpoint - 0.1, y_midpoint + 0.1], color="black", lw=1)
        ax.plot([x1, x1], [y_midpoint - 0.1, y_midpoint + 0.1], color="black", lw=1)

        scale_distance_human_readable = get_human_readable_number_of_bp(scale_distance)

        ax.text(
            (x0 + (scale_distance / 2)),
            y_midpoint - 0.2,
            scale_distance_human_readable,
            ha="center",
            va="center",
        )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)


class GenomicAxis(cb.XAxis):
    """
    A genomic axis that shows the genomic coordinates of the region.
    """

    def __init__(self, **kwargs):
        super(GenomicAxis, self).__init__()
        self.properties.update(kwargs)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax

        ax.set_xlim(gr.start, gr.end)
        ticks = np.linspace(gr.start, gr.end, 10)
        labels = [f"{x:,.0f}" for x in ticks]

        ax.axis["x"] = ax.new_floating_axis(0, 0.5)
        ax.axis["x"].axis.set_ticks(ticks)
        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis["x"].axis.set_tick_params(which="minor", bottom="on")

        ax.axis["x"].major_ticklabels.set(size=10)

        if "where" in self.properties and self.properties["where"] == "top":
            ax.axis["x"].set_axis_direction("top")


class MatrixCapcruncherAverage(MatrixCapcruncher):
    """
    A matrix track that averages multiple matrices.

    Args:
        matricies (List[MatrixCapcruncher]): List of MatrixCapcruncher objects
        aggregation (str, optional): Aggregation method. Defaults to "mean". Choices are "sum", "mean", "median".
        **kwargs: Additional arguments to pass to the track

    """

    def __init__(
        self,
        matricies: List[MatrixCapcruncher],
        aggregation: Literal["sum", "mean", "median"] = "mean",
        **kwargs,
    ):
        self.matricies = matricies
        self.aggregation = aggregation
        self.properties = matricies[0].properties
        self.properties.update(kwargs)
        self.properties["name"] = f"CCMatrix.{self.properties.get('title')}"

        # Need to override the coolbox default if we need a cmap to be set
        self.properties["color"] = kwargs.get("color", self.properties["color"])

        # Override the defaults
        self.properties["balance"] = "no"

    @property
    def subtracks(self):
        return self.matricies

    @functools.cache
    def fetch_data(self, gr: cb.GenomeRange, gr2=None, **kwargs):
        data = []
        for matrix in tqdm.tqdm(self.matricies):
            data.append(matrix.fetch_data(gr, **kwargs))

        try:
            func_agg = getattr(np, self.aggregation)
        except AttributeError:
            raise ValueError(
                f"Aggregation function {self.aggregation} not found in numpy"
            )

        # Aggregate the list of matricies into a single matrix
        data = func_agg(data, axis=0)

        self.fetched_gr = gr
        self.fetched_gr2 = gr2
        return data

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        # fetch processed plot_data
        self.matrix = self.fetch_data(gr, **kwargs)
        # plot matrix
        img = self.plot_matrix(gr, kwargs.get("gr2"))
        self.adjust_figure(gr, kwargs.get("gr2"))
        self.draw_colorbar(img)
        self.plot_label()


class HighlightsFromFile(cb.HighLightsFromFile):
    """
    Plot highlights from a bed file.
    This is a modified version of the HighlightsFromFile class that does not print block of highlights but instead highlights the entire region.

    Args:
        file (os.PathLike): Path to the bed file
        **kwargs: Additional arguments to pass to the track
    """

    def plot(self, ax, gr: cb.GenomeRange, **kwargs):
        from matplotlib.patches import Rectangle

        regions = self.fetch_data(gr, **kwargs)

        for start, end, color in regions:
            if self.properties["color"] != "bed_rgb":
                color = self.properties["color"]

            ymin, ymax = ax.get_ylim()

            # Add a small offset to the y-axis to avoid the highlights overlapping with the axis
            # ymin += 0.001 * (ymax - ymin)
            # ymax -= 0.001 * (ymax - ymin)

            ax.add_patch(
                Rectangle(
                    (start, ymin),
                    end - start,
                    ymax - ymin,
                    color=color,
                    alpha=self.properties["alpha"],
                    linewidth=0,
                )
            )
