import os
import pathlib
from enum import Enum
from typing import Any, Dict, List, Literal, Optional, Union

import coolbox.api as cb
import numpy as np

from plotnado.api.tracks import (
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


class TrackType(Enum):
    heatmap = MatrixCapcruncher
    heatmap_summary = MatrixCapcruncherAverage
    bigwig_fragment = BigwigFragment
    bigwig_fragment_collection = BigwigFragmentCollection
    bigwig_fragment_collection_overlay = BigwigFragmentCollectionOverlay
    bigwig = cb.BigWig
    bigwig_summary = BigwigOverlay
    bed_memory = BedMemory
    bed_simple = BedSimple
    xaxis = GenomicAxis
    scale = ScaleBar


class FilelessTracks(Enum):
    spacer = cb.Spacer
    scale = ScaleBar
    xaxis = GenomicAxis


class TrackWrapper:
    """
    Provides a wrapper around tracks to provide a consistent interface

    Args:
        file (os.PathLike): Path to file to plot
        track_class: Class to use for the track or alias
        **kwargs: Additional arguments to pass to the track
    """

    def __init__(
        self,
        track_type: Union[str, TrackType],
        file: Optional[Union[str, List[str]]] = None,
        **kwargs,
    ):
        """
        Initialize a TrackWrapper

        Args:
            track_type (Union[str, TrackType]): Type of track to plot
            file (Optional[Union[str, List[str]]], optional): Path to file to plot. Defaults to None.
            **kwargs: Additional arguments to pass to the track
        """

        self.track_type = track_type
        self.file = file
        self.properties = dict()
        self.properties.update(kwargs)

    @property
    def track_class(self):
        """
        Get the track class
        """

        try:
            track_class = TrackType[self.track_type].value
        except KeyError:
            if getattr(cb, self.track_type):
                track_class = getattr(cb, self.track_type)
            else:
                raise ValueError(
                    f"Unknown track type {self.track_type}, select from: {', '.join([t.name for t in TrackType])}"
                )

        return track_class

    def get_track(self) -> cb.Track:
        """
        Get the track object with the specified properties and adding track type specific properties
        """

        track = self.track_class(self.file, **self.properties)
        return track

    @property
    def path(self) -> str:
        if isinstance(self.file, (list, tuple, np.ndarray)):
            return [str(pathlib.Path(f).resolve()) for f in self.file]
        else:
            return str(pathlib.Path(self.file).resolve())

    def __repr__(self) -> str:
        return f"WrappedTrack({self.properties.get('title')}, {self.properties.get('type')})"


class Figure:
    """
    Generates a figure from a list of tracks

    Args:
        tracks (List[TrackWrapper], optional): List of tracks to plot. Defaults to None.
        auto_spacing (bool, optional): Automatically add a spacer track between each track. Defaults to False.
        **kwargs: Additional arguments to pass to the figure
    """

    def __init__(
        self,
        tracks: List[TrackWrapper] = None,
        auto_spacing: bool = False,
        frame_args: Dict[str, Any] = None,
        **kwargs,
    ) -> None:

        self.frame = cb.Frame(**frame_args)
        self.auto_spacing = auto_spacing
        self.properties = dict()
        self.properties.update(kwargs)

        if tracks:
            self.tracks = set(tracks)
            self.add_tracks(tracks)
        else:
            self.tracks = set()

    def add_track(self, track: TrackWrapper) -> None:
        """
        Add a track to the figure

        Args:
            track (TrackWrapper): Track to add
        """
        self.tracks.add(track)
        self.frame.add_track(track.get_track())

    def add_tracks(self, tracks: List[TrackWrapper]) -> None:
        """
        Add a list of tracks to the figure

        Args:
            tracks (List[TrackWrapper]): List of tracks to add
        """
        for track in tracks:
            if self.auto_spacing:
                spacer = TrackWrapper(None, track_type="spacer")
                self.add_track(spacer.get_track())

            self.add_track(track)

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

        tracks = []
        for track_name, attr in config.items():
            file = attr.pop("file") if attr.get("file") else None
            track_name = attr.pop("title") if attr.get("title") else track_name
            tracks.append(TrackWrapper(file, title=track_name, **attr))
        return cls(tracks, **kwargs)

    @classmethod
    def from_frame(cls, frame: cb.Frame, **kwargs) -> "Figure":
        """
        Instantiate a Figure from a coolbox Frame

        Args:
            frame (cb.Frame): coolbox Frame to instantiate from
            **kwargs: Additional arguments to pass to the figure
        """
        tracks = []
        for track in frame.tracks:
            tracks.append(TrackWrapper(track.properties["file"], **track.properties))

        return cls(tracks, **kwargs)

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

        def _get_n_tracks_of_type(config: Dict[str, Dict], track_type: str):
            return sum(1 for t in config.keys() if track_type in t)

        config = OrderedDict()
        for track in self.tracks:
            # Perform conversions for file-less tracks
            if track.properties.get("type") in ["spacer", "scale", "xaxis"]:
                track_type = track.properties.get("type")
                n = _get_n_tracks_of_type(config, track_type)
                config[f"{track_type} {n}"] = track.properties
                config[f"{track_type} {n}"]["file"] = None
            elif track.properties.get("type") == "genes":
                track_type = track.properties.get("type")
                n = _get_n_tracks_of_type(config, track_type)
                config[f"{track_type} {n}"] = track.properties
                config[f"{track_type} {n}"]["file"] = track.path
            else:
                config[track.properties["title"]] = track.properties
                config[track.properties["title"]]["file"] = track.path

        outfile = output if output else "config.toml"

        with open(outfile, "w") as f:
            config_str = toml.dumps(config)
            f.write(config_str)

        if not output:
            return config
