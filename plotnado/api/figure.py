import os
import pathlib
from collections import OrderedDict, defaultdict
from enum import Enum
from typing import Any, Dict, List, Literal, Optional, Union

import matplotlib.patches

import coolbox.api as cb
import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import matplotlib.figure

from .track_wrapper import FilelessTracks, TrackType, TrackWrapper
from .tracks import (
    Autoscaler,
    BedMemory,
    BedSimple,
    BigwigFragment,
    BigwigFragmentCollection,
    BigwigFragmentCollectionOverlay,
    BigwigOverlay,
    GenomicAxis,
    MatrixCapcruncher,
    MatrixCapcruncherAverage,
    ScaleBar,
    HighlightsFromFile,
)


def get_track_title(track: TrackWrapper, config: Dict[str, Dict]) -> str:
    """
    Get the title of a track

    Args:
        track (TrackWrapper): Track to get the title of
        config (Dict[str, Dict]): Configuration of the tracks

    Returns:
        str: Title of the track
    """

    try:
        title = track.properties.get("title", track.track_type)
    except AttributeError:
        title = ""

    if title in config:
        n = 1
        while f"{title} {n}" in config:  # Ensure the title is unique
            n += 1
        title = f"{title} {n}"

    return title


def plot_label_on_track(
    track,
    ax,
    ymin,
    ymax,
    gr: cb.GenomeRange,
    location: Literal["left", "right"] = "left",
):
    ydelta = ymax - ymin
    small_x = 0.01 * gr.length

    if not hasattr(track, 'collection'):
        if location == "left":
            ax.text(
                gr.start,
                ymax - ydelta * 0.2,
                track.properties["title"],
                horizontalalignment="left",
                verticalalignment="top",
            )
        elif location == "right":
            ax.text(
                gr.end - small_x,
                ymax - ydelta * 0.2,
                track.properties["title"],
                horizontalalignment="right",
                verticalalignment="top",
            )
    
    else:

        y = ymax - ydelta * 0.2
        y_rec_top = y + (y * 0.1)

        if location == "left":
            for i, bw in enumerate(track.collection):

                rec_start = gr.start
                rec_end = rec_start + (gr.length * 0.001)
                tiny_x = 0.0001 * gr.length

                ax.add_patch(
                    matplotlib.patches.Polygon(
                        [
                            (rec_start, y - (i * 0.1)),
                            (rec_start, y_rec_top - (i * 0.1)),
                            (rec_end, y_rec_top - (i * 0.1)),
                            (rec_end , y - (i * 0.1)),
                        ],
                        color=bw.properties.get("color", "black"),
                    )
                )

                track.label_ax.text(
                    rec_end + tiny_x,
                    0.55 - (i * 0.1),
                    bw.properties.get("title", f"Subtrack {i}"),
                    horizontalalignment="left",
                    verticalalignment="center",
                )
        
        elif location == "right":
            for i, bw in enumerate(track.collection):

                rec_start = gr.end - (gr.length * 0.001)
                rec_end = rec_start + (gr.length * 0.001)
                tiny_x = 0.0001 * gr.length

                ax.add_patch(
                    matplotlib.patches.Polygon(
                        [
                            (rec_start, y - (i * 0.1)),
                            (rec_start, y_rec_top - (i * 0.1)),
                            (rec_end, y_rec_top - (i * 0.1)),
                            (rec_end , y - (i * 0.1)),
                        ],
                        color=bw.properties.get("color", "black"),
                    )
                )

                track.label_ax.text(
                    rec_start - tiny_x,
                    0.55 - (i * 0.1),
                    bw.properties.get("title", f"Subtrack {i}"),
                    horizontalalignment="right",
                    verticalalignment="center",
                )


def _plot(self, *args, label_loc: Literal["left", "right"] = "right", **kwargs):
    """
    Plot all tracks.

    >>> from coolbox.api import *
    >>> frame = XAxis() + XAxis()
    >>> frame.plot("chr1", 100000, 200000)
    >>> frame.plot("chr1:100000-200000")
    """
    from copy import copy
    from collections import OrderedDict

    import matplotlib
    import matplotlib.pyplot as plt
    import mpl_toolkits.axisartist as axisartist

    from coolbox.utilities import (
        cm2inch,
        GenomeRange,
        get_logger,
    )

    log = get_logger(__name__)

    if len(args) == 3:
        gr = GenomeRange(*args[:3])
        gr2 = None
    elif len(args) == 1:
        gr = GenomeRange(args[0])
        gr2 = None
    elif len(args) == 6:
        gr = GenomeRange(args[:3])
        gr2 = GenomeRange(args[3:])
    elif len(args) == 2:
        gr = GenomeRange(args[0])
        gr2 = GenomeRange(args[1])
    elif self.current_range:
        gr = self.current_range
        gr2 = None
    else:
        raise ValueError(
            "Please specify a genomic range in uscs format. For example: 'chr1:100000-200000'"
        )
    # cache for the previous GenomeRange
    self.goto(gr, gr2)

    tracks_height = self.get_tracks_height()
    self.properties["height"] = sum(tracks_height)

    fig = plt.figure(
        figsize=cm2inch(self.properties["width"], self.properties["height"])
    )
    if "title" in self.properties:
        fig.suptitle(self.properties["title"])

    grids = matplotlib.gridspec.GridSpec(
        len(tracks_height),
        3,
        height_ratios=tracks_height,
        width_ratios=self.properties["width_ratios"],
        wspace=0.01,
    )

    axis_list = []
    for idx, track in enumerate(self.tracks.values()):
        if label_loc == "left":
            # Label -> Y-axis -> Track
            label_ax = plt.subplot(grids[idx, 0])
            label_ax.set_axis_off()

            y_ax = plt.subplot(grids[idx, 1])
            y_ax.set_axis_off()

            ax = axisartist.Subplot(fig, grids[idx, 2])
            fig.add_subplot(ax)
            ax.axis[:].set_visible(False)
            ax.patch.set_visible(False)

        else:
            # Y-axis -> Track -> Label
            y_ax = plt.subplot(grids[idx, 0])
            y_ax.set_axis_off()

            ax = axisartist.Subplot(fig, grids[idx, 1])
            fig.add_subplot(ax)
            ax.axis[:].set_visible(False)
            ax.patch.set_visible(False)

            label_ax = plt.subplot(grids[idx, 2])
            label_ax.set_axis_off()

        track.label_ax = label_ax
        track.y_ax = y_ax

        try:
            # Attention, copy is necessary, otherwise GenomeRange may change due to call of gr.change_chrom_names
            track.plot(ax, copy(gr), gr2=copy(gr2))

            if track.properties.get("label_on_track") in ['True', 'yes', 'T', 'Y', '1', True, 1]:
                ymin, ymax = track.adjust_plot(ax, gr)
                plot_label_on_track(track, ax, ymin, ymax, gr, location=track.properties['label_loc'])

        except Exception as e:
            import sys, os

            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = exc_tb.tb_frame.f_code.co_filename
            log.error(
                "Error occured when plot track:\n"
                "\ttrack name: {}\n\ttrack type:{}\n"
                "\tError: {} {}\n"
                '\toccurred in "{}", line {}'.format(
                    track.name, type(track), type(e), str(e), fname, exc_tb.tb_lineno
                )
            )
            log.exception(e)
        track.plot_coverages(ax, gr, gr2)

        axis_list.append(ax)

    margins = self.properties["margins"]
    fig.subplots_adjust(
        wspace=0,
        hspace=0.0,
        left=margins["left"],
        right=margins["right"],
        bottom=margins["bottom"],
        top=margins["top"],
    )

    plt.close()

    return fig


class Figure:
    """
    Generates a figure from a list of tracks

    Args:
        tracks (List[TrackWrapper], optional): List of tracks to plot. Defaults to None.
        auto_spacing (bool, optional): Automatically add a spacer track between each track. Defaults to False.
        frame_args (Dict[str, Any], optional): Arguments to pass to the frame. Defaults to None.
        autoscale_groups (Dict[str, List[int]], optional): Groups of tracks to autoscale together. Defaults to None.
        **kwargs: Additional arguments to pass to the figure
    """

    def __init__(
        self,
        tracks: List[TrackWrapper] = None,
        autospacing: bool = False,
        autospacing_height: float = 0.1,
        frame_args: Dict[str, Any] = None,
        highlight_regions: Union[pathlib.Path, cb.HighLights] = None,
        highlight_regions_color: str = "blue",
        highlight_regions_kwargs: Dict[str, Any] = None,
        autocolor: bool = False,
        **kwargs,
    ) -> None:
        self.frame = cb.Frame(**frame_args if frame_args else dict())
        self.autospacing = autospacing
        self.autospacing_height = autospacing_height
        self.autocolor = autocolor
        self.properties = dict()
        self.properties.update(kwargs)

        self.tracks = OrderedDict()
        if tracks:
            self.add_tracks(tracks)

        if isinstance(highlight_regions, (pathlib.Path, str)):
            self.highlight_regions = HighlightsFromFile(
                highlight_regions,
                color=highlight_regions_color,
                **highlight_regions_kwargs,
            )
        else:
            self.highlight_regions = highlight_regions

        self.track_name_to_frame_name_mapping = dict()

    def add_track(self, track: Union[str, TrackWrapper], **kwargs) -> None:
        """
        Add a track to the figure

        Args:
            track (TrackWrapper): Track to add
        """

        if isinstance(track, str):
            track = TrackWrapper(track_type=track, **kwargs)

        # Add a spacer track if auto_spacing is enabled
        if self.autospacing:
            spacer = TrackWrapper(track_type="spacer", height=self.autospacing_height)
            title = get_track_title(spacer, self.tracks)
            self.tracks[title] = spacer

        # Get the title of the track
        title = get_track_title(track, self.tracks)

        # Add the track to the collection
        self.tracks[title] = track

        # Add the track to the frame so coolbox can plot its
        self.frame.add_track(track.track)

        # Add the track to the mapping (take the last track added to the frame as the frame name)
        self.track_name_to_frame_name_mapping[title] = list(self.frame.tracks.keys())[
            -1
        ]

    def add_tracks(self, tracks: List[TrackWrapper]) -> None:
        """
        Add a list of tracks to the figure

        Args:
            tracks (List[TrackWrapper]): List of tracks to add
        """
        for track in tracks:
            self.add_track(track)

    def _autoscale(self, gr: cb.GenomeRange, gr2: cb.GenomeRange = None):
        # Extract the autoscale groups
        scale_groups = defaultdict(list)
        for title, track in self.tracks.items():
            if track.autoscale_group:
                scale_groups[track.autoscale_group].append(title)

        # Autoscale the tracks
        for group, tracks in scale_groups.items():
            # Need to translate the track names to the frame names
            tracks = [self.track_name_to_frame_name_mapping[title] for title in tracks]
            autoscaler = Autoscaler(
                [self.frame.tracks[title] for title in tracks],
                gr,
                gr2,
            )

            max_value = autoscaler.max_value
            min_value = autoscaler.min_value

            for title in tracks:
                self.frame.tracks[title].properties["max_value"] = max_value
                self.frame.tracks[title].properties["min_value"] = min_value

    def _autocolor(self):
        colors = list(plt.cm.tab20.colors)
        tracktypes_for_autocolor = ["bigwig"]

        for track_name, track in self.frame.tracks.items():
            if any(y in track_name.lower() for y in tracktypes_for_autocolor):
                color_number = np.random.choice(range(len(colors)))
                color = colors[color_number]
                track.properties["color"] = color

    @property
    def data_tracks(self):
        tracks = OrderedDict()
        for title, track in self.tracks.items():
            if not any(y in title for y in ["spacer"]):
                tracks[title] = track
        return tracks

    def plot(
        self,
        gr: Union[str, cb.GenomeRange],
        gr2: Union[str, cb.GenomeRange] = None,
        show: bool = True,
        extend: Union[int, Dict[Literal["upstream", "downstream"], int]] = 0,
        **kwargs,
    ) -> matplotlib.figure.Figure:
        """
        Plot the figure

        Args:
            gr (Union[str, GenomeRange]): GenomeRange to plot
            gr2 (Union[str, GenomeRange], optional): Second GenomeRange to plot. Defaults to None.
            show (bool, optional): Show the figure. Defaults to True.
            **kwargs: Additional arguments to pass to the plot

        Returns:
            matplotlib.figure.Figure: The figure
        """

        # Ensure the genome ranges are valid
        if isinstance(gr, str):
            gr = cb.GenomeRange(gr)

        if isinstance(gr2, str):
            gr2 = cb.GenomeRange(gr2)

        # Extend the genome ranges
        if isinstance(extend, int):
            extend = {"upstream": extend, "downstream": extend}

        gr.start = max(0, gr.start - extend["upstream"])
        gr.end = gr.end + extend["downstream"]

        if gr2:
            gr2.start = max(0, gr2.start - extend["upstream"])
            gr2.end = gr2.end + extend["downstream"]

        # Autoscale the tracks
        self._autoscale(gr, gr2)

        # Autocolor the tracks if specified
        if self.autocolor:
            self._autocolor()

        # Highlight the regions if specified
        if self.highlight_regions:
            self.frame = self.frame * self.highlight_regions

        # Patch the frame plot method to allow for label location
        import functools

        self.frame.plot = functools.partial(_plot, self.frame)

        # Plot the figure as 2D or 1D depending on the number of genome ranges
        if gr2:
            fig = self.frame.plot(gr, gr2, **kwargs)
        else:
            fig = self.frame.plot(gr, **kwargs)
        if show:
            fig.show()

        return fig

    def plot_gene(
        self,
        gene: str,
        genome: str,
        **kwargs,
    ):
        """
        Plot the figure for a gene

        Args:
            gene (str): Gene to plot
            genome (str): Genome to plot the gene on
            **kwargs: Additional arguments to pass to the plot
        """
        import importlib
        import json

        genes_prefix = importlib.resources.files("plotnado.data.genes")
        with open(genes_prefix / f"{genome}.json") as f:
            genes = json.load(f)

        gene = genes[gene]
        gr = cb.GenomeRange(gene["chrom"], gene["start"], gene["end"])
        return self.plot(gr, **kwargs)

    def plot_regions(
        self,
        regions: Union[
            pathlib.Path, pr.PyRanges, pd.DataFrame, Dict[str, cb.GenomeRange]
        ],
        show: bool = True,
        **kwargs,
    ) -> Dict[str, matplotlib.figure.Figure]:
        """
        Plot the figure for a list of regions

        Args:
            regions (Union[str, pd.DataFrame, List[GenomeRange]]): Regions to plot
            show (bool, optional): Show the figure. Defaults to True.
            **kwargs: Additional arguments to pass to the plot

        Returns:
            Dict[str, matplotlib.figure.Figure]: Dictionary of figures for each region
        """

        # Format the regions into a dictionary of genome ranges
        if isinstance(regions, pathlib.Path):
            regions = pr.read_bed(str(regions))

        if isinstance(regions, pr.PyRanges):
            regions = {
                (
                    row.Name if row.Name else f"{row.Chromosome}_{row.Start}_{row.End}"
                ): cb.GenomeRange(
                    chrom=row.Chromosome,
                    start=row.Start,
                    end=row.End,
                )
                for _, row in regions.df.itertuples()
            }

        if isinstance(regions, pd.DataFrame):
            assert (
                "Chromosome" in regions.columns
            ), "Chromosome column not found in regions"
            assert "Start" in regions.columns, "Start column not found in regions"
            assert "End" in regions.columns, "End column not found in regions"

            regions = {
                (
                    row.Name if row.Name else f"{row.Chromosome}_{row.Start}_{row.End}"
                ): cb.GenomeRange(
                    chrom=row.Chromosome,
                    start=row.Start,
                    end=row.End,
                )
                for _, row in regions.iterrows()
            }

        plots = dict()
        for name, gr in regions.items():
            plots[name] = self.plot(gr, show=False, **kwargs)

        return plots

    def save(
        self,
        gr: Union[str, cb.GenomeRange],
        gr2: Union[str, cb.GenomeRange] = None,
        output: str = None,
        **kwargs,
    ) -> None:
        """
        Plots the figure and saves it to a file

        Args:
            gr (Union[str, GenomeRange]): GenomeRange to plot
            gr2 (Union[str, GenomeRange], optional): Second GenomeRange to plot. Defaults to None.
            output (str, optional): Path to save the figure to. Defaults to None.
            **kwargs: Additional arguments to pass to the plot
        """

        fig = self.plot(gr, gr2, show=False, **kwargs)
        if output:
            fig.savefig(output, dpi=300)
        else:
            fig.savefig(f"{gr.chrom}_{gr.start}_{gr.end}.png", dpi=300)

    @classmethod
    def from_toml(cls, toml_file: os.PathLike, **kwargs) -> "Figure":
        """
        Instantiate a Figure from a toml file

        Args:
            toml_file (os.PathLike): Path to toml file
            **kwargs: Additional arguments to pass to the figure
        """
        import toml

        with open(toml_file) as f:
            config = toml.load(f)

        _cls = cls()
        for track_name, track_properties in config.items():
            _cls.add_track(TrackWrapper.from_dict(track_properties))

        return _cls

    def to_toml(self, output: str = None) -> Union[None, Dict[str, Any]]:
        """
        Save the Figure to a toml file

        Args:
            output (str, optional): Path to save the toml file to. Defaults to None.

        Returns:
            Union[None, Dict[str, Any]]: If output is not specified, returns a dict of the toml file

        """

        from collections import OrderedDict

        import toml

        # Ordered dict with the key being a unique identifier for the track
        # ideally the title of the track or the track type with a number to ensure uniqueness
        config = OrderedDict()
        for title, track in self.tracks.items():
            config[title] = track.to_dict()

        outfile = output if output else "config.toml"

        with open(outfile, "w") as f:
            config_str = toml.dumps(config)
            f.write(config_str)

        if not output:
            return config

    def __repr__(self) -> str:
        return f"Figure({len(self.tracks)} tracks)"
