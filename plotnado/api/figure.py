import os
import pathlib
from collections import OrderedDict
from enum import Enum
from typing import Any, Dict, List, Literal, Optional, Union

import coolbox.api as cb
import numpy as np
import pandas as pd

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

    title = track.properties.get("title", track.track_type)
    if title in config:
        n = 1
        while f"{title} {n}" in config: # Ensure the title is unique
            n += 1
        title = f"{title} {n}"

    return title


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
        auto_spacing: bool = False,
        frame_args: Dict[str, Any] = None,
        autoscale_groups: Dict[str, List[int]] = None,
        **kwargs,
    ) -> None:

        self.frame = cb.Frame(**frame_args if frame_args else dict())
        self.auto_spacing = auto_spacing
        self.autoscale_groups = autoscale_groups
        self.properties = dict()
        self.properties.update(kwargs)

        self.tracks = OrderedDict()
        if tracks:
            self.add_tracks(tracks)

    def add_track(self, track: TrackWrapper) -> None:
        """
        Add a track to the figure

        Args:
            track (TrackWrapper): Track to add
        """
        # Add a spacer track if auto_spacing is enabled
        if self.auto_spacing:
            spacer = TrackWrapper(track_type="spacer")
            title = get_track_title(spacer, self.tracks)
            self.tracks[title] = spacer

        # Get the title of the track
        title = get_track_title(track, self.tracks)

        # Add the track to the collection
        self.tracks[title] = track

        # Add the track to the frame so coolbox can plot its
        self.frame.add_track(track.track)

    def add_tracks(self, tracks: List[TrackWrapper]) -> None:
        """
        Add a list of tracks to the figure

        Args:
            tracks (List[TrackWrapper]): List of tracks to add
        """
        for track in tracks:
            self.add_track(track)

    def _autoscale(self, gr: cb.GenomeRange, gr2: cb.GenomeRange = None):
        # Deal with autoscaling
        if self.autoscale_groups:
            for group, tracks_indexes in self.autoscale_groups.items():
                tracks_for_scaling = []
                for index, track in enumerate(self.frame.tracks.values()):
                    if index in tracks_indexes:
                        tracks_for_scaling.append(track)

                autoscaler = Autoscaler(tracks_for_scaling, gr, gr2)
                for track in tracks_for_scaling:
                    track.properties["max_value"] = autoscaler.max_value
                    track.properties["min_value"] = autoscaler.min_value

    def plot(
        self,
        gr: Union[str, cb.GenomeRange],
        gr2: Union[str, cb.GenomeRange] = None,
        show: bool = True,
        **kwargs,
    ) -> None:
        """
        Plot the figure

        Args:
            gr (Union[str, GenomeRange]): GenomeRange to plot
            gr2 (Union[str, GenomeRange], optional): Second GenomeRange to plot. Defaults to None.
            show (bool, optional): Show the figure. Defaults to True.
            **kwargs: Additional arguments to pass to the plot
        """

        self._autoscale(gr, gr2)

        if gr2:
            fig = self.frame.plot(gr, gr2, **kwargs)
        else:
            fig = self.frame.plot(gr, **kwargs)
        if show:
            fig.show()

        return fig

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
            print(config)
            config_str = toml.dumps(config)
            f.write(config_str)

        if not output:
            return config

    def __repr__(self) -> str:
        return f"Figure({len(self.tracks)} tracks)"
