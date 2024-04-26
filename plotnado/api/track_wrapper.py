import pathlib
from collections import OrderedDict
from enum import Enum
import re
from typing import Any, Dict, List, Optional, Tuple, Union

import coolbox.api as cb
import numpy as np

from . import (
    BedMemory,
    BedSimple,
    BigwigFragment,
    BigwigFragmentCollection,
    BigwigFragmentCollectionOverlay,
    BigwigOverlay,
    Genes,
    GenomicAxis,
    MatrixCapcruncher,
    MatrixCapcruncherAverage,
    ScaleBar,
    BigwigSubtraction
)


URL_REGEX = re.compile(
    r"http[s]?://"  # http:// or https://
    r"(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|"  # domain...
    r"localhost|"  # localhost...
    r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
    r"(?::\d+)?"  # optional port
    r"(?:/?|[/?]\S+)$",
    re.IGNORECASE,
)

MATRIX_TRACKS = (MatrixCapcruncher, MatrixCapcruncherAverage, cb.Cool)
BIGWIG_TRACKS = (
    BigwigFragment,
    BigwigFragmentCollection,
    BigwigFragmentCollectionOverlay,
    BigwigOverlay,
    BigwigSubtraction,
    cb.BigWig,
)
CUSTOM_TRACKS = (BedMemory, BedSimple, GenomicAxis, ScaleBar)

AGGREGATED_TRACKS_MEMMORY = (
    MatrixCapcruncherAverage,
    BigwigOverlay,
    BigwigFragmentCollection,
)
AGGREGATED_TRACKS_FILE = (BigwigFragmentCollectionOverlay, BigwigSubtraction)

ALLOWED_NON_FILE_TRACKS = (
    GenomicAxis,
    ScaleBar,
    BedMemory,
    BedSimple,
    Genes,
    cb.HighLightsFromFile,
)

HAS_PRESETS = (Genes)


class TrackType(Enum):
    heatmap = MatrixCapcruncher
    heatmap_summary = MatrixCapcruncherAverage
    bigwig_fragment = BigwigFragment
    bigwig_fragment_collection = BigwigFragmentCollection
    bigwig_fragment_collection_overlay = BigwigFragmentCollectionOverlay
    bigwig_subtraction = BigwigSubtraction
    bigwig = cb.BigWig
    bigwig_overlay = BigwigOverlay
    bed_memory = BedMemory
    bed_simple = BedSimple
    xaxis = GenomicAxis
    scale = ScaleBar
    scalebar = ScaleBar
    spacer = cb.Spacer
    genes = Genes


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
        track_type: Union[str, TrackType, Any],
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

        # Check if the track type is a string, enum or class
        # and convert it to a string

        # First check if the track type is an enum
        if isinstance(track_type, TrackType):
            self.track_type = track_type.name

        # If it is a string, check if it is a valid track type
        elif isinstance(track_type, str):
            # Check if the track type is a valid track type
            if track_type in [t.name for t in TrackType]:
                self.track_type = TrackType[track_type].name
            elif hasattr(cb, track_type):
                # Is the track type a class in the coolbox api
                self.track_type = track_type

        elif hasattr(cb, type(track_type).__name__):
            # Check if the track type is a class in the coolbox api
            self.track_type = type(track_type).__name__

        else:
            # If the track type is not a valid track type
            # raise an error
            raise ValueError(
                f"Unknown track type {track_type}, select from: {', '.join([t.name for t in TrackType])} or provide a custom track class"
            )

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
            elif self.track_type in [*CUSTOM_TRACKS, *BIGWIG_TRACKS, *MATRIX_TRACKS]:
                track_class = self.track_type
            else:
                raise ValueError(
                    f"Unknown track type {self.track_type}, select from: {', '.join([t.name for t in TrackType])} or provide a custom track class"
                )

        return track_class

    @property
    def track(self) -> cb.Track:
        """
        Get the track object with the specified properties and adding track type specific properties
        """

        if self.file is None and self.track_class in [t.value for t in FilelessTracks]:
            track = self.track_class(**self.properties)

        elif self.track_class in AGGREGATED_TRACKS_FILE:
            assert all(
                [pathlib.Path(p).exists() for p in self.file]
            ), f"Not files provided to track {self.track_class} exist"
            track = self.track_class(self.file, **self.properties)

        elif self.track_class in AGGREGATED_TRACKS_MEMMORY:
            tracks = []
            for t in self.file:
                if isinstance(t, TrackWrapper):
                    tracks.append(t.track)
                else:
                    tracks.append(t)

            track = self.track_class(tracks, **self.properties)

        elif self.track_class in ALLOWED_NON_FILE_TRACKS:
            track = self.track_class(self.file, **self.properties)

        elif pathlib.Path(self.file).exists():
            track = self.track_class(self.file, **self.properties)

        elif self.properties.get("ignore_file_validation"):
            track = self.track_class(self.file, **self.properties)
        
        elif URL_REGEX.match(self.file):
            track = self.track_class(self.file, **self.properties)

        else:
            raise FileNotFoundError(f"File {self.file} does not exist")

        return track

    @property
    def path(self) -> str:
        """
        Get the path to the 'file' attribute these can be either a ssingle file or a list of files
        """
        _files = []
        if isinstance(self.file, (list, tuple, np.ndarray)):

            for f in self.file:
                if isinstance(f, (TrackWrapper, cb.Track)):
                    if hasattr(f, "file"):
                        _file = f.file
                    elif f.properties.get("file"):
                        _file = f.properties["file"]
                elif isinstance(f, (str, pathlib.Path)):
                    _file = f
                _files.append(_file)
        elif self.file is None:
            return None
        else:
            _files.append(self.file)
        
        files = []
        for f in _files:
            if URL_REGEX.match(f):
                files.append(f)
            elif pathlib.Path(f).exists():
                files.append(str(pathlib.Path(f).resolve()))
            elif self.properties.get("ignore_file_validation") or self.track.properties.get("ignore_file_validation"):
                files.append(f)
            else:
                raise FileNotFoundError(f"File {f} does not exist")
        
        if len(files) == 1:
            return files[0]
        else:
            return tuple(files)
    
    @property
    def autoscale_group(self) -> str:
        """
        Get the autoscale group
        """
        return self.properties.get("autoscale_group", None)
            

    def to_dict(self) -> Dict:
        """
        Convert the TrackWrapper to a dictionarys
        """

        # Check to see if the track has subtracks
        # need to handle this differently
        if hasattr(self.track, "subtracks"):
            subtracks = OrderedDict()
            for track in self.track.subtracks:
                subtracks[track.properties["name"]] = TrackWrapper(
                    track_type=type(track).__name__, **track.properties
                ).to_dict()

            return {
                "track_type": self.track_type,
                "file": self.path,
                "subtracks": subtracks,
                **self.properties,
            }
        else:
            return {
                "track_type": self.track_type,
                "file": self.path,
                **self.properties,
            }

    @classmethod
    def from_dict(cls, data: Dict):
        """
        Load the TrackWrapper from a dictionary
        """

        if "subtracks" in data:
            subtracks = []
            for subtrack_name, subtrack_data in data["subtracks"].items():
                subtracks.append(TrackWrapper.from_dict(subtrack_data).track)
            data["file"] = subtracks

        return cls(**data)

    def __repr__(self) -> str:
        title = self.properties.get("title", "No title")
        track_class = self.track_class
        return f"TrackWrapper({title}, {track_class})"
