"""Figure class for composing and plotting genomic tracks."""

from __future__ import annotations

import importlib.resources
import json
import math
from pathlib import Path
from typing import Self

import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger

from .tracks import (
    BedTrack,
    BigWigCollection,
    BigWigDiff,
    BigWigTrack,
    CapcruncherTrack,
    CoolerAverage,
    CoolerTrack,
    Genes,
    GenomicAxis,
    GenomicRegion,
    HighlightsFromFile,
    HLineTrack,
    LinksTrack,
    NarrowPeakTrack,
    ScaleBar,
    Spacer,
    Track,
    VLineTrack,
    BigwigOverlay,
)


class Figure:
    """Compose and plot multiple genomic tracks."""

    def __init__(
        self,
        tracks: list[Track] | None = None,
        width: float = 12,
        track_height: float = 2.0,
    ):
        self.tracks: list[Track] = tracks or []
        self.width = width
        self.track_height = track_height
        self._highlight_regions: list[GenomicRegion] = []
        self._autocolor_palette: str | None = None
        self._autoscale: bool = False

    def add_track(self, track: str | Track, **kwargs) -> Self:
        if isinstance(track, str):
            track = self._create_track_from_alias(track, **kwargs)
        self.tracks.append(track)
        return self

    def autoscale(self, enable: bool = True) -> Self:
        self._autoscale = enable
        return self

    def autocolor(self, palette: str = "tab10") -> Self:
        self._autocolor_palette = palette
        import matplotlib.colors as mcolors

        cmap = plt.get_cmap(palette)
        for index, track in enumerate(self.tracks):
            if hasattr(track, "aesthetics") and hasattr(track.aesthetics, "color"):
                track.aesthetics.color = mcolors.to_hex(cmap(index % cmap.N))
        return self

    def highlight(self, region: str | GenomicRegion) -> Self:
        self._highlight_regions.append(GenomicRegion.into(region))
        return self

    def _create_track_from_alias(self, alias: str, **kwargs) -> Track:
        alias_map = {
            "scalebar": ScaleBar,
            "scale": ScaleBar,
            "genes": Genes,
            "spacer": Spacer,
            "bigwig": BigWigTrack,
            "bed": BedTrack,
            "axis": GenomicAxis,
            "highlight": HighlightsFromFile,
            "bigwig_overlay": BigwigOverlay,
            "narrowpeak": NarrowPeakTrack,
            "links": LinksTrack,
            "hline": HLineTrack,
            "vline": VLineTrack,
            "cooler": CoolerTrack,
            "capcruncher": CapcruncherTrack,
            "cooler_average": CoolerAverage,
            "bigwig_collection": BigWigCollection,
            "bigwig_diff": BigWigDiff,
        }
        key = alias.lower()
        if key not in alias_map:
            raise ValueError(
                f"Unknown track alias: {alias}. Available: {list(alias_map.keys())}"
            )
        return alias_map[key](**kwargs)

    @staticmethod
    def _apply_extend(gr: GenomicRegion, extend: float | int | None) -> GenomicRegion:
        if not extend:
            return gr

        if isinstance(extend, float):
            if extend < 0:
                raise ValueError("extend as float must be >= 0")
            bp = int(gr.length * extend)
        else:
            bp = int(extend)

        if bp < 0:
            raise ValueError("extend must be >= 0")
        return gr.extend(upstream=bp, downstream=bp)

    def _apply_autoscale_groups(
        self, axes: list[matplotlib.axes.Axes], main_tracks: list[Track]
    ) -> None:
        grouped: dict[str, list[int]] = {}
        for index, track in enumerate(main_tracks):
            if track.autoscale_group:
                grouped.setdefault(track.autoscale_group, []).append(index)

        for indices in grouped.values():
            mins = []
            maxs = []
            for index in indices:
                y_limits = axes[index].get_ylim()
                mins.append(y_limits[0])
                maxs.append(y_limits[1])
            group_min = min(mins)
            group_max = max(maxs)
            if group_min == group_max:
                group_max = group_min + 1
            for index in indices:
                axes[index].set_ylim(group_min, group_max)

    def plot(
        self,
        region: str | GenomicRegion,
        show: bool = True,
        extend: float | int | None = None,
        **kwargs,
    ) -> matplotlib.figure.Figure | None:
        gr = GenomicRegion.into(region)
        gr = self._apply_extend(gr, extend)

        if not self.tracks:
            logger.warning("No tracks to plot")
            return None

        main_tracks = [track for track in self.tracks if track.height > 0]
        global_tracks = [track for track in self.tracks if track.height == 0]
        heights = [track.height for track in main_tracks]
        total_height = max(1.0, sum(heights) * self.track_height)

        fig = plt.figure(figsize=(self.width, total_height), **kwargs)
        fig.patch.set_facecolor("white")

        from matplotlib.gridspec import GridSpec

        gs = GridSpec(len(main_tracks), 1, height_ratios=heights, hspace=0.1)
        axes = [fig.add_subplot(gs[index]) for index in range(len(main_tracks))]

        def create_overlay_ax(zorder: int, label: str):
            ax = fig.add_axes(axes[0].get_position(), label=label, zorder=zorder)
            pos_top = axes[0].get_position()
            pos_bottom = axes[-1].get_position()
            ax.set_position(
                [
                    pos_top.x0,
                    pos_bottom.y0,
                    pos_top.width,
                    pos_top.y1 - pos_bottom.y0,
                ]
            )
            ax.patch.set_alpha(0)
            from .tracks.utils import clean_axis

            clean_axis(ax)
            ax.set_xlim(gr.start, gr.end)
            return ax

        bg_ax = create_overlay_ax(zorder=-1, label="background_overlay")
        fg_ax = create_overlay_ax(zorder=10, label="foreground_overlay")

        if self._autoscale:
            from .tracks.scaling import Autoscaler

            Autoscaler(tracks=self.tracks, gr=gr).apply()

        for ax, track in zip(axes, main_tracks):
            ax.patch.set_alpha(0)
            try:
                track.plot(ax, gr)
            except Exception as exc:
                logger.error(f"Error plotting track {track.__class__.__name__}: {exc}")
                ax.text(0.5, 0.5, f"Error: {exc}", ha="center", va="center", transform=ax.transAxes)

        self._apply_autoscale_groups(axes, main_tracks)

        for highlight_gr in self._highlight_regions:
            if highlight_gr.chromosome == gr.chromosome:
                start = max(highlight_gr.start, gr.start)
                end = min(highlight_gr.end, gr.end)
                if start < end:
                    bg_ax.axvspan(start, end, color="#ffd700", alpha=0.15, zorder=-1)

        for track in global_tracks:
            try:
                if isinstance(track, HighlightsFromFile):
                    track.plot(bg_ax, gr)
                elif isinstance(track, VLineTrack):
                    track.plot(fg_ax, gr)
                elif hasattr(track, "plot_on_axes"):
                    track.plot_on_axes(gr, axes)
            except Exception as exc:
                logger.error(f"Error plotting global track {track.__class__.__name__}: {exc}")

        fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)
        if show:
            plt.show()
        return fig

    def _resolve_gene_track(self) -> Genes:
        for track in self.tracks:
            if isinstance(track, Genes):
                return track
        return Genes(genome="hg38")

    def plot_gene(self, gene: str, extend: float = 0.5, **kwargs) -> matplotlib.figure.Figure | None:
        genes_track = self._resolve_gene_track()
        has_genes_track = any(isinstance(track, Genes) for track in self.tracks)

        if genes_track.genome:
            bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
            with open(bed_prefix / "genes.json") as handle:
                mapping = json.load(handle)
            gene_file = bed_prefix / mapping[genes_track.genome]
            genes_df = pd.read_csv(gene_file, sep="\t", header=None)
        elif genes_track.data:
            genes_df = pd.read_csv(str(genes_track.data), sep="\t", header=None)
        else:
            raise ValueError("No gene annotation source available for plot_gene")

        genes_df.columns = [
            "chrom",
            "start",
            "end",
            "name",
            *[f"field_{index}" for index in range(max(0, genes_df.shape[1] - 4))],
        ]
        match = genes_df.loc[genes_df["name"].astype(str).str.upper() == gene.upper()]
        if match.empty:
            raise ValueError(f"Gene {gene} not found in annotation source")

        row = match.iloc[0]
        region = GenomicRegion(chromosome=row["chrom"], start=int(row["start"]), end=int(row["end"]))

        if has_genes_track:
            return self.plot(region, extend=extend, **kwargs)

        # If no genes track is present, include one for this plotting call so
        # the gene context is visible by default.
        self.tracks.append(genes_track)
        try:
            return self.plot(region, extend=extend, **kwargs)
        finally:
            self.tracks.pop()

    def plot_regions(
        self, regions: list[str] | str, ncols: int = 1, **kwargs
    ) -> list[matplotlib.figure.Figure]:
        if ncols < 1:
            raise ValueError("ncols must be >= 1")

        show = kwargs.pop("show", True)
        region_strings: list[str]
        path_candidate = Path(regions) if isinstance(regions, str) else None

        if isinstance(regions, str) and path_candidate.exists():
            bed_df = pd.read_csv(path_candidate, sep="\t", header=None, comment="#")
            region_strings = [
                f"{row[0]}:{int(row[1])}-{int(row[2])}" for _, row in bed_df.iterrows()
            ]
        elif isinstance(regions, str):
            region_strings = [regions]
        else:
            region_strings = regions

        figs: list[matplotlib.figure.Figure] = []
        for region in region_strings:
            fig = self.plot(region, show=False, **kwargs)
            if fig is not None:
                figs.append(fig)

        if not figs:
            return []

        if ncols > 1 and len(figs) > 1:
            nrows = math.ceil(len(figs) / ncols)
            grid_fig, grid_axes = plt.subplots(
                nrows,
                ncols,
                figsize=(self.width * ncols, self.track_height * max(2, nrows * 2)),
            )
            axes_array = np.atleast_1d(grid_axes).ravel()

            for index, source_fig in enumerate(figs):
                source_fig.canvas.draw()
                image = np.asarray(source_fig.canvas.buffer_rgba())
                axes_array[index].imshow(image)
                axes_array[index].set_axis_off()
                axes_array[index].set_title(region_strings[index], fontsize=9)

            for index in range(len(figs), len(axes_array)):
                axes_array[index].set_axis_off()

            for source_fig in figs:
                plt.close(source_fig)

            figs = [grid_fig]

        if show:
            for fig in figs:
                fig.show()
        return figs

    @staticmethod
    def _track_registry() -> dict[str, type[Track]]:
        return {
            cls.__name__: cls
            for cls in [
                BedTrack,
                BigWigTrack,
                BigWigCollection,
                BigWigDiff,
                BigwigOverlay,
                Genes,
                GenomicAxis,
                HighlightsFromFile,
                HLineTrack,
                LinksTrack,
                NarrowPeakTrack,
                ScaleBar,
                Spacer,
                VLineTrack,
                CoolerTrack,
                CapcruncherTrack,
                CoolerAverage,
            ]
        }

    def to_toml(self, path: str) -> None:
        tracks_by_type: dict[str, list[dict]] = {}
        for track in self.tracks:
            tracks_by_type.setdefault(track.__class__.__name__, []).append(track.model_dump())

        payload = {
            "figure": {"width": self.width, "track_height": self.track_height},
            "tracks": tracks_by_type,
        }

        try:
            import tomli_w
        except ImportError as exc:
            raise ImportError("to_toml requires optional dependency 'tomli-w'") from exc

        with open(path, "wb") as handle:
            handle.write(tomli_w.dumps(payload).encode("utf-8"))

    @classmethod
    def from_toml(cls, path: str) -> "Figure":
        try:
            import tomllib
        except ImportError:
            import tomli as tomllib  # type: ignore[no-redef]

        with open(path, "rb") as handle:
            payload = tomllib.load(handle)

        figure_data = payload.get("figure", {})
        fig = cls(
            width=figure_data.get("width", 12),
            track_height=figure_data.get("track_height", 2.0),
        )

        registry = cls._track_registry()
        tracks_payload = payload.get("tracks", [])

        if isinstance(tracks_payload, list):
            for track_spec in tracks_payload:
                type_name = track_spec["type"]
                params = track_spec.get("params", {})
                if type_name not in registry:
                    raise ValueError(f"Unknown track type in TOML: {type_name}")
                fig.add_track(registry[type_name](**params))
            return fig

        if isinstance(tracks_payload, dict):
            for type_name, raw_params in tracks_payload.items():
                if type_name not in registry:
                    raise ValueError(f"Unknown track type in TOML: {type_name}")

                params_list = raw_params if isinstance(raw_params, list) else [raw_params]
                for params in params_list:
                    fig.add_track(registry[type_name](**(params or {})))
            return fig

        raise ValueError("Invalid TOML format: 'tracks' must be a list or table")

    def save(
        self,
        path: str | Path,
        region: str | GenomicRegion,
        dpi: int = 300,
        **kwargs,
    ) -> None:
        fig = self.plot(region, show=False)
        if fig:
            fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
            logger.info(f"Saved figure to {path}")
            plt.close(fig)

    def __repr__(self) -> str:
        track_names = [track.__class__.__name__ for track in self.tracks]
        return f"Figure(tracks={track_names})"
