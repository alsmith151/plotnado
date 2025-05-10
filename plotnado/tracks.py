from abc import ABC, abstractmethod
from collections import OrderedDict, defaultdict
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import matplotlib.axes
import matplotlib.axis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandera as pa
import pybigtools
import seaborn as sns
from loguru import logger
from pydantic import BaseModel, Field
from pandera.typing import Index, DataFrame, Series
import matplotlib


def clean_axis(ax: matplotlib.axes.Axes) -> None:
    """
    Clean the axis:
    - Remove the ticks
    - Remove the spines

    Args:
        ax (matplotlib.axes.Axes): The axis object
    """
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
    ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())





class GenomicRegion(BaseModel):
    chromosome: str
    start: int
    end: int
    strand: Literal["+", "-"] = "+"


    @property
    def length(self) -> int:
        return self.end - self.start
    
    @property
    def center(self) -> int:
        return (self.start + self.end) // 2


    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}({self.strand})"

    @classmethod
    def from_str(cls, region_str: str) -> "GenomicRegion":
        chromosome, region = region_str.split(":")
        start_end, strand = region.split("(")
        start, end = start_end.split("-")
        return cls(
            chromosome=chromosome, start=int(start), end=int(end), strand=strand[:-1]
        )

    @classmethod
    def from_tuple(cls, region_tuple: tuple) -> "GenomicRegion":
        return cls(
            chromosome=region_tuple[0],
            start=region_tuple[1],
            end=region_tuple[2],
            strand=region_tuple[3],
        )

    @classmethod
    def from_list(cls, region_list: list) -> "GenomicRegion":
        return cls(
            chromosome=region_list[0],
            start=region_list[1],
            end=region_list[2],
            strand=region_list[3],
        )

    @classmethod
    def from_dict(cls, region_dict: dict) -> "GenomicRegion":
        return cls(
            chromosome=region_dict["chromosome"],
            start=region_dict["start"],
            end=region_dict["end"],
            strand=region_dict["strand"],
        )

    @classmethod
    def from_named_tuple(cls, region_named_tuple: tuple) -> "GenomicRegion":
        return cls(
            chromosome=region_named_tuple.chromosome,
            start=region_named_tuple.start,
            end=region_named_tuple.end,
            strand=region_named_tuple.strand,
        )

    @classmethod
    def into(cls, gr: Union[str , tuple , list , dict]) -> "GenomicRegion":
        try:
            match type(gr):
                case "str":
                    return cls.from_str(gr)
                case "tuple":
                    return cls.from_tuple(gr)
                case "list":
                    return cls.from_list(gr)
                case "dict":
                    return cls.from_dict(gr)
                case "namedtuple":
                    return cls.from_named_tuple(gr)
                case _:
                    raise ValueError("Unsupported type")
        except ValueError as e:
            logger.error(e)
            raise e


class BigwigAesthetics(BaseModel):

    # Plot style
    style: Literal["std", "scatter", "line", "heatmap"] = "std"

    # General aesthetics
    color: str = "black"
    fill: bool = True
    alpha: float = 1.0

    # Scatter plot aesthetics
    scatter_point_size: float = 1.0

    # Hist plot aesthetics
    linewidth: float = 1.0

    min_value: Optional[float] = None
    max_value: Optional[float] = None

    plot_title: bool = True
    title_location: Literal["left", "right"] = "left"
    title_height: float = 0.5 # fraction of the y-axis
    
    plot_scale: bool = True
    scale_location: Literal["left", "right"] = "left"
    scale_height: float = 0.5 # fraction of the y-axis



class Aesthetics(Enum):
    BIGWIG: BigwigAesthetics


class Track(BaseModel, ABC):
    title: str
    data: Path | pd.DataFrame
    aesthetics: Aesthetics

    class Config:
        arbitrary_types_allowed = True

    def fetch_data(self) -> pd.DataFrame:
        pass

    def plot(self) -> None:
        pass

    def save(self, path: Path) -> None:
        pass



class Labeller(BaseModel):
    gr: GenomicRegion
    y_min: float
    y_max: float

    plot_title: bool = True
    plot_scale: bool = True

    title: str
    title_size: int = 8
    title_color: str = "black"
    title_font: str = "Arial"
    title_weight: Literal["normal", "bold"] = "bold"
    title_location: Literal["left", "right"] = "left"
    title_height: float = 0.8 # fraction of the y-axis

    scale_min: float
    scale_max: float
    scale_precision: int = 2 # number of decimal places
    scale_size: int = 6
    scale_color: str = "black"
    scale_font: str = "Arial"
    scale_weight: Literal["normal", "bold"] = "bold"
    scale_location: Literal["left", "right"] = "right"
    scale_height: float = 0.8 # fraction of the y-axis


    @property
    def y_delta(self) -> float:
        return self.y_max - self.y_min

    def _plot_title(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """
        Plot the title of the track

        Args:
            ax (matplotlib.axes.Axes): The axis object
            gr (GenomicRegion): The genomic region
        """

        if self.title_location == "left":
            ax.text(
                gr.start + (0.01 * gr.length),
                self.y_delta * self.title_height,
                self.title,
                horizontalalignment="left",
                verticalalignment="top",
                fontdict=
                {
                    "size": self.title_size,
                    "color": self.title_color,
                    "fontname": self.title_font,
                    "weight": self.title_weight,
                }
            )
        elif self.title_location == "right":
            ax.text(
                gr.end - (0.01 * gr.length),
                self.y_delta * self.title_height,
                self.title,
                horizontalalignment="right",
                verticalalignment="top",
                fontdict=
                {
                    "size": self.title_size,
                    "color": self.title_color,
                    "fontname": self.title_font,
                    "weight": self.title_weight,
                }
            )

    
    def _format_scale(self, value: float) -> str:
        if value % 1 == 0:
            return str(int(value))
        else:
            return f"{value:.{self.scale_precision}f}"''
    
    def _plot_scale(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:

        y_min = self._format_scale(self.y_min)
        y_max = self._format_scale(self.y_max)


        if self.scale_location == "right":
            ax.text(
                gr.end - (0.01 * gr.length),
                self.y_delta * self.scale_height,
                f"[ {y_min} - {y_max} ]",
                horizontalalignment="right",
                verticalalignment="top",
                fontdict={
                    "size": self.scale_size,
                    "color": self.scale_color,
                    "fontname": self.scale_font,
                    "weight": self.scale_weight,
                }
            )
        elif self.scale_location == "left":
            ax.text(
                gr.start + (0.01 * gr.length),
                self.y_delta * self.scale_height,
                f"[ {y_min} - {y_max} ]",
                horizontalalignment="left",
                verticalalignment="top",
                fontdict={
                    "size": self.scale_size,
                    "color": self.scale_color,
                    "fontname": self.scale_font,
                    "weight": self.scale_weight,
                }
            )

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        if self.plot_title:
            self._plot_title(ax, gr)
        if self.plot_scale:
            self._plot_scale(ax, gr)
        
        # Clean the axis
        clean_axis(ax)

        



class BedgraphDataFrame(pa.DataFrameModel):
    chrom: Series[str]
    start: Series[int]
    end: Series[int]
    value: Series[float]

    class Config:
        coerce = True
        strict = True


class BigWigTrack(Track):
    aesthetics: BigwigAesthetics

    y_min: Optional[float] = None
    y_max: Optional[float] = None

    def _fetch_from_disk(self, gr: GenomicRegion) -> BedgraphDataFrame:
        bw = pybigtools.open(self.data)
        records = bw.records(gr.chromosome, gr.start, gr.end)
        df = pd.DataFrame(records, columns=["start", "end", "value"]).assign(
            chrom=gr.chromosome
        )
        return BedgraphDataFrame(df)

    def _fetch_from_df(self, gr: GenomicRegion) -> pd.DataFrame:
        return BedgraphDataFrame(
            self.data[
                (self.data["chrom"] == gr.chromosome)
                & (self.data["start"] >= gr.start)
                & (self.data["end"] <= gr.end)
            ]
        )

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:
        if isinstance(self.data, Path):
            data = self._fetch_from_disk(gr)
        elif isinstance(self.data, pd.DataFrame):
            data = self._fetch_from_df(gr)

        return data

    def _plot_stairs(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, values: BedgraphDataFrame
    ) -> None:
        edges = np.linspace(gr.start, gr.end, values.shape[0] + 1)
        ax.stairs(
            edges=edges,
            values=values["value"],
            linewidth=self.aesthetics.linewidth,
            color=self.aesthetics.color,
            alpha=self.aesthetics.alpha,
            fill=self.aesthetics.fill,
        )
    

    def _plot_scatter(
        self, ax: matplotlib.axes.Axes, gr: GenomicRegion, values: BedgraphDataFrame
    ) -> None:
        

        x = np.linspace(gr.start, gr.end, values.shape[0])
        ax.scatter(
            x,
            values["value"],
            color=self.aesthetics.color,
            alpha=self.aesthetics.alpha,
            s=self.aesthetics.scatter_point_size,
        )



    def plot(self, gr: GenomicRegion, ax: matplotlib.axes.Axes) -> None:
        
        # Fetch the data
        data = self.fetch_data(gr)

        # Plot the data
        if self.aesthetics.style == 'std':
            self._plot_stairs(ax, gr, data)
        
        elif self.aesthetics.style == 'scatter':
            self._plot_scatter(ax, gr, data)
        

        # Fix the limits of the plot
        ax.set_xlim(gr.start, gr.end)
        self.y_min = data['value'].min() if self.aesthetics.min_value is None else self.aesthetics.min_value
        self.y_max = data['value'].max() if self.aesthetics.max_value is None else self.aesthetics.max_value 
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(ymin=self.y_min, ymax=self.y_max)

        # Plot the labels
        labeller = Labeller(
            gr=gr,
            y_min=self.y_min,
            y_max=self.y_max,
            plot_title=self.aesthetics.plot_title,
            plot_scale=self.aesthetics.plot_scale,
            title=self.title,
            scale_min=self.y_min,
            scale_max=self.y_max
        )

        labeller.plot(ax, gr)








