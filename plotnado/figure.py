"""GenomicFigure class for composing and plotting genomic tracks."""

from __future__ import annotations

import importlib.resources
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Self

import matplotlib.axes
import matplotlib.figure
import matplotlib.lines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger
from pydantic import BaseModel

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
    OverlayTrack,
    ScaleBar,
    Spacer,
    Track,
    VLineTrack,
    BigwigOverlay,
    LabelConfig,

    QuantNadoCoverageTrack,
    QuantNadoStrandedCoverageTrack,
    QuantNadoMethylationTrack,
    QuantNadoVariantTrack,
    list_options,
)
from .theme import BuiltinTheme, Theme
from .tracks.registry import registry
from .figure_methods import GenomicFigureMethods

if TYPE_CHECKING:
    from .template import Template


class GenomicFigure(GenomicFigureMethods):
    """Compose and plot multiple genomic tracks.

    Example:
        fig = (
            GenomicFigure()
            .bigwig("signal.bw", title="H3K27ac")
            .genes()
            .axis()
            .scalebar()
        )
        fig.plot("chr1:1000000-2000000")
    """

    def __init__(
        self,
        tracks: list[Track] | None = None,
        width: float = 12,
        track_height: float = 1.25,
        theme: Theme | BuiltinTheme | str | None = BuiltinTheme.PUBLICATION,
    ):
        """Create a figure container for genomic tracks.

        Args:
            tracks: Optional pre-populated list of tracks.
            width: Figure width in inches.
            track_height: Height scaling factor per track.
            theme: Optional `Theme` (or builtin name like "default", "minimal",
                "publication") applied to tracks on add/initialization.
        """
        self.tracks: list[Track] = tracks or []
        self.width = width
        self.track_height = track_height
        self.theme = self._resolve_theme(theme)
        self.highlight_color = (
            self.theme.highlight_color if self.theme is not None else "#ffd700"
        )
        self.highlight_alpha = self.theme.highlight_alpha if self.theme is not None else 0.15
        self._highlight_regions: list[GenomicRegion] = []
        self._autocolor_palette: str | list[str] | None = None
        self._autoscale: bool = False

        if self.theme is not None:
            for track in self.tracks:
                self.theme.apply(track)
            self._apply_theme_palette()

    @classmethod
    def from_template(
        cls,
        template: Template | str | Path,
        *,
        width: float | None = None,
        theme: Theme | BuiltinTheme | str | None = BuiltinTheme.PUBLICATION,
    ) -> GenomicFigure:
        """Build a GenomicFigure from a template file or Template object.

        Args:
            template: A ``Template`` instance, or a path to a YAML template file.
            width: Override figure width in inches. Defaults to the template's width.
            theme: Theme to apply. Defaults to the publication theme.

        Returns:
            A fully configured GenomicFigure ready for ``.plot()`` or ``.save()``.

        Example::

            fig = GenomicFigure.from_template("template.yaml")
            fig.save("out.png", region="chr1:1000000-2000000")
        """
        from .template import Template as _Template
        from .render import TemplateCompiler

        if not isinstance(template, _Template):
            template = _Template.load(template)

        plan = TemplateCompiler.compile(template)
        fig = cls(
            width=width if width is not None else plan.width,
            track_height=plan.track_height,
            theme=theme,
        )

        if plan.add_scalebar:
            fig.scalebar()
        if plan.add_axis:
            fig.axis()
        if plan.add_genes and plan.genome:
            fig.genes(plan.genome)

        for resolved in plan.tracks:
            kwargs = resolved.to_figure_kwargs()
            data = resolved.get_data()
            track_type_str = str(resolved.track_spec.type)
            if data is not None:
                fig.add_track(track_type_str, **{resolved.source_kwarg_name(): data}, **kwargs)
            else:
                fig.add_track(track_type_str, **kwargs)

        return fig

    @staticmethod
    def _resolve_theme(theme: Theme | BuiltinTheme | str | None) -> Theme | None:
        if theme is None:
            return None
        if isinstance(theme, Theme):
            return theme
        try:
            return Theme.from_builtin(theme)
        except ValueError as exc:
            raise ValueError(
                f"Unknown builtin theme: {theme}. Available: {[item.value for item in BuiltinTheme]}"
            ) from exc

    def add_track(self, track: str | Track, **kwargs: Any) -> Self:
        """Add a track instance or track alias to the figure.

        Args:
            track: Track instance or alias (for example `"bigwig"`, `"genes"`).
            **kwargs: Parameters used when `track` is provided as an alias.

        Returns:
            Self, enabling method chaining.
        """
        if isinstance(track, str):
            track = self._create_track_from_alias(track, **kwargs)

        # Apply theme defaults to the track
        if self.theme is not None:
            self.theme.apply(track)

        self.tracks.append(track)
        if self._autocolor_palette is not None:
            self._apply_autocolor()
        else:
            self._apply_theme_palette()
        return self

    @staticmethod
    def _aesthetic_explicitly_set(track: Track, field_name: str) -> bool:
        aesthetics = getattr(track, "aesthetics", None)
        explicit_fields = getattr(aesthetics, "model_fields_set", set())
        return field_name in explicit_fields

    @staticmethod
    def _set_auto_color(track: Track, color: str) -> None:
        track.color = color
        aesthetics = getattr(track, "aesthetics", None)
        explicit_fields = getattr(aesthetics, "model_fields_set", None)
        if isinstance(explicit_fields, set):
            explicit_fields.discard("color")

    @staticmethod
    def _label_field_explicitly_set(track: Track, field_name: str) -> bool:
        label = getattr(track, "label", None)
        explicit_fields = getattr(label, "model_fields_set", set())
        return field_name in explicit_fields

    @staticmethod
    def _is_theme_palette_eligible(track: Track) -> bool:
        if not track.has_aesthetic("color"):
            return False
        if isinstance(
            track,
            (
                ScaleBar,
                GenomicAxis,
                HighlightsFromFile,
                Spacer,
                HLineTrack,
                VLineTrack,
                Genes,
                CoolerTrack,
                CapcruncherTrack,
                CoolerAverage,
            ),
        ):
            return False
        return True

    def _apply_theme_palette(self) -> None:
        if self.theme is None or not self.theme.auto_palette or not self.theme.palette:
            return

        eligible_tracks = [track for track in self.tracks if self._is_theme_palette_eligible(track)]
        if len(eligible_tracks) < 2:
            return

        auto_color_tracks = [
            track
            for track in eligible_tracks
            if self._is_theme_palette_eligible(track) and not self._aesthetic_explicitly_set(track, "color")
        ]
        if not auto_color_tracks:
            return

        explicit_group_colors = {
            track.color_group: track.color
            for track in eligible_tracks
            if track.color_group and self._aesthetic_explicitly_set(track, "color")
        }

        assigned_groups: dict[str, str] = {}
        palette_index = 0
        for track in auto_color_tracks:
            group = track.color_group
            if group and group in explicit_group_colors:
                self._set_auto_color(track, explicit_group_colors[group])
                continue
            if group and group in assigned_groups:
                self._set_auto_color(track, assigned_groups[group])
                continue

            color = self.theme.palette[palette_index % len(self.theme.palette)]
            palette_index += 1
            self._set_auto_color(track, color)
            if group:
                assigned_groups[group] = color

    @staticmethod
    def _is_meta_track(track: Track) -> bool:
        return isinstance(
            track,
            (
                ScaleBar,
                GenomicAxis,
                HighlightsFromFile,
                Spacer,
                HLineTrack,
                VLineTrack,
            ),
        )

    @classmethod
    def _should_autocolor_track(cls, track: Track) -> bool:
        return track.has_aesthetic("color") and not cls._is_meta_track(track)

    @classmethod
    def available_track_aliases(cls) -> dict[str, str]:
        """Return available alias -> TrackClass mappings."""
        return {
            alias: entry.cls.__name__
            for alias, entry in registry.all_entries().items()
        }

    @classmethod
    def track_options(cls, track: str | type[Track]) -> dict[str, dict]:
        """Return programmatic option metadata for a track alias or class.

        Examples:
            GenomicFigure.track_options("bigwig")
            GenomicFigure.track_options(BigWigTrack)
        """
        if isinstance(track, str):
            key = track.lower()
            try:
                entry = registry.get(key)
                track_cls = entry.cls
            except KeyError:
                available = sorted(registry.all_entries().keys())
                raise ValueError(
                    f"Unknown track alias: {track}. Available: {available}"
                ) from None
        else:
            track_cls = track
        return list_options(track_cls)

    @classmethod
    def track_options_markdown(cls, track: str | type[Track]) -> str:
        """Return markdown-formatted option tables for a track alias or class."""
        if isinstance(track, str):
            key = track.lower()
            try:
                entry = registry.get(key)
                track_cls = entry.cls
            except KeyError:
                available = sorted(registry.all_entries().keys())
                raise ValueError(
                    f"Unknown track alias: {track}. Available: {available}"
                ) from None
        else:
            track_cls = track
        return track_cls.options_markdown()

    def autoscale(self, enable: bool = True) -> Self:
        """Enable or disable automatic y-axis autoscaling across tracks."""
        self._autoscale = enable
        return self

    def autocolor(self, palette: str | list[str] | None = None) -> Self:
        """Apply a matplotlib palette across tracks that expose a `color` field."""
        if palette is None:
            if self.theme is not None and self.theme.palette:
                palette = list(self.theme.palette)
            else:
                palette = "tab10"

        self._autocolor_palette = palette
        self._apply_autocolor()
        return self

    def _apply_autocolor(self) -> None:
        if self._autocolor_palette is None:
            return

        import matplotlib.colors as mcolors

        auto_tracks = [
            track
            for track in self.tracks
            if self._should_autocolor_track(track)
        ]
        if not auto_tracks:
            return

        explicit_group_colors = {
            track.color_group: track.color
            for track in auto_tracks
            if track.color_group and self._aesthetic_explicitly_set(track, "color")
        }
        assigned_groups: dict[str, str] = {}
        color_index = 0

        for track in auto_tracks:
            group = track.color_group
            if group and group in explicit_group_colors:
                self._set_auto_color(track, explicit_group_colors[group])
                continue
            if group and group in assigned_groups:
                self._set_auto_color(track, assigned_groups[group])
                continue

            if isinstance(self._autocolor_palette, list):
                if not self._autocolor_palette:
                    continue
                color = self._autocolor_palette[color_index % len(self._autocolor_palette)]
            else:
                cmap = plt.get_cmap(self._autocolor_palette)
                color = mcolors.to_hex(cmap(color_index % cmap.N))
            color_index += 1
            self._set_auto_color(track, color)
            if group:
                assigned_groups[group] = color

    def highlight(self, region: str | GenomicRegion) -> Self:
        """Register a genomic region to draw as a background highlight."""
        self._highlight_regions.append(GenomicRegion.into(region))
        return self

    def highlight_style(self, color: str | None = None, alpha: float | None = None) -> Self:
        """Set default color/alpha used by highlighted regions."""
        if color is not None:
            self.highlight_color = color
        if alpha is not None:
            self.highlight_alpha = alpha
        return self

    def _create_track_from_alias(self, alias: str, **kwargs: Any) -> Track:
        key = alias.lower()
        try:
            entry = registry.get(key)
        except KeyError:
            available = sorted([e.track_type.value for e in registry.all_entries().values()])
            raise ValueError(
                f"Unknown track type: {alias!r}. Available: {available}"
            ) from None
        track_cls = entry.cls

        track_fields = set(track_cls.model_fields.keys())
        aesthetics_model = track_cls.aesthetics_model()
        aesthetics_fields = (
            set(aesthetics_model.model_fields.keys()) if aesthetics_model is not None else set()
        )
        label_fields = set(LabelConfig.model_fields.keys())

        track_kwargs: dict[str, Any] = {}
        aesthetics_overrides: dict[str, Any] = {}
        label_overrides: dict[str, Any] = {}

        for field_name, field_value in kwargs.items():
            if field_name in {"aesthetics", "label"} or field_name in track_fields:
                track_kwargs[field_name] = field_value
            elif field_name in aesthetics_fields:
                aesthetics_overrides[field_name] = field_value
            elif field_name in label_fields:
                label_overrides[field_name] = field_value
            else:
                track_kwargs[field_name] = field_value

        if aesthetics_model is not None and aesthetics_overrides:
            base_aesthetics = track_kwargs.get("aesthetics")
            if isinstance(base_aesthetics, BaseModel):
                payload = base_aesthetics.model_dump()
            elif isinstance(base_aesthetics, dict):
                payload = dict(base_aesthetics)
            elif base_aesthetics is None:
                payload = {}
            else:
                payload = {}
            payload.update(aesthetics_overrides)
            track_kwargs["aesthetics"] = aesthetics_model(**payload)

        if label_overrides:
            base_label = track_kwargs.get("label")
            if isinstance(base_label, BaseModel):
                label_payload = base_label.model_dump()
            elif isinstance(base_label, dict):
                label_payload = dict(base_label)
            elif base_label is None:
                label_payload = {}
            else:
                label_payload = {}
            label_payload.update(label_overrides)
            track_kwargs["label"] = LabelConfig(**label_payload)

        return track_cls(**track_kwargs)

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

    def _font_family(self) -> str:
        if self.theme is not None:
            if self.theme.font_family:
                return self.theme.font_family
            return self.theme.label.title_font
        return "DejaVu Sans"

    def plot(
        self,
        region: str | GenomicRegion,
        show: bool = False,
        extend: float | int | None = None,
        **kwargs: Any,
    ) -> matplotlib.figure.Figure | None:
        """Render all tracks for a single genomic region.

        Args:
            region: Region string (`chr:start-end`) or `GenomicRegion`.
            show: Whether to call `plt.show()` before returning.
            extend: Optional symmetric extension (fraction or base pairs).
            **kwargs: Forwarded to `matplotlib.figure.Figure` creation.

        Returns:
            The rendered matplotlib figure, or `None` if no tracks are present.
        """
        gr = GenomicRegion.into(region)
        gr = self._apply_extend(gr, extend)
        # Ensure every track receives the exact same plotting bounds.
        gr = GenomicRegion(chromosome=gr.chromosome, start=int(gr.start), end=int(gr.end))

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

        hspace = self.theme.subplot_hspace if self.theme is not None else 0.08
        gs = GridSpec(len(main_tracks), 1, height_ratios=heights, hspace=hspace)
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
                ax.set_xlim(gr.start, gr.end)
            except Exception as exc:
                logger.error(f"Error plotting track {track.__class__.__name__}: {exc}")
                ax.text(0.5, 0.5, f"Error: {exc}", ha="center", va="center", transform=ax.transAxes)

        self._apply_autoscale_groups(axes, main_tracks)

        for highlight_gr in self._highlight_regions:
            if highlight_gr.chromosome == gr.chromosome:
                start = max(highlight_gr.start, gr.start)
                end = min(highlight_gr.end, gr.end)
                if start < end:
                    bg_ax.axvspan(
                        start,
                        end,
                        color=self.highlight_color,
                        alpha=self.highlight_alpha,
                        zorder=-1,
                    )

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

        font_family = self._font_family()
        for axis in fig.axes:
            for text_artist in axis.texts:
                text_artist.set_fontfamily(font_family)

        if self.theme is not None:
            fig.subplots_adjust(
                left=self.theme.margin_left,
                right=self.theme.margin_right,
                top=self.theme.margin_top,
                bottom=self.theme.margin_bottom,
            )
            if self.theme.separator_color and len(axes) > 1:
                left = axes[0].get_position().x0
                right = axes[0].get_position().x1
                for index, axis in enumerate(axes[:-1]):
                    track_above = main_tracks[index]
                    track_below = main_tracks[index + 1]
                    if isinstance(track_above, Genes):
                        continue
                    if (
                        isinstance(track_above, ScaleBar)
                        and isinstance(track_below, Genes)
                    ) or (
                        isinstance(track_above, Genes)
                        and isinstance(track_below, ScaleBar)
                    ):
                        continue
                    y = axis.get_position().y0
                    separator = matplotlib.lines.Line2D(
                        [left, right],
                        [y, y],
                        transform=fig.transFigure,
                        color=self.theme.separator_color,
                        alpha=self.theme.separator_alpha,
                        linewidth=self.theme.separator_linewidth,
                        zorder=20,
                    )
                    fig.add_artist(separator)
        else:
            fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)

        if show:
            plt.show()
        else:
            # Detach from pyplot state to avoid duplicate auto-rendering in notebooks
            # when the returned Figure object is also displayed as cell output.
            plt.close(fig)
        return fig

    def _resolve_gene_track(self) -> Genes:
        for track in self.tracks:
            if isinstance(track, Genes):
                return track
        return Genes(genome="hg38")

    def plot_gene(self, gene: str, extend: float = 0.5, **kwargs: Any) -> matplotlib.figure.Figure | None:
        """Resolve a gene symbol to a region and plot it.

        Args:
            gene: Gene symbol to look up in configured gene annotations.
            extend: Fractional extension around gene bounds.
            **kwargs: Forwarded to `plot`.

        Returns:
            The rendered figure or `None` when plotting is skipped.
        """
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
        self, regions: list[str] | str, ncols: int = 1, **kwargs: Any
    ) -> list[matplotlib.figure.Figure]:
        """Plot one or many regions, optionally composing a multi-column grid.

        Args:
            regions: A region string, list of region strings, or BED-like path.
            ncols: Number of columns for grid composition when >1 region.
            **kwargs: Forwarded to `plot`.

        Returns:
            A list of matplotlib figures.
        """
        if ncols < 1:
            raise ValueError("ncols must be >= 1")

        show = kwargs.pop("show", True)
        region_strings: list[str]
        path_candidate = Path(regions) if isinstance(regions, str) else None

        if isinstance(regions, str) and path_candidate is not None and path_candidate.exists():
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
                OverlayTrack,
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

                QuantNadoCoverageTrack,
                QuantNadoStrandedCoverageTrack,
                QuantNadoMethylationTrack,
                QuantNadoVariantTrack,
            ]
        }

    def to_toml(self, path: str) -> None:
        """Serialize figure and tracks to a TOML file."""
        def _prune_none(value):
            if isinstance(value, dict):
                return {
                    key: _prune_none(item)
                    for key, item in value.items()
                    if item is not None
                }
            if isinstance(value, list):
                return [_prune_none(item) for item in value if item is not None]
            return value

        tracks_by_type: dict[str, list[dict]] = {}
        for track in self.tracks:
            tracks_by_type.setdefault(track.__class__.__name__, []).append(
                _prune_none(track.model_dump())
            )

        payload = {
            "figure": {"width": self.width, "track_height": self.track_height},
            "tracks": tracks_by_type,
        }
        payload = _prune_none(payload)

        try:
            import tomli_w
        except ImportError as exc:
            raise ImportError("to_toml requires optional dependency 'tomli-w'") from exc

        with open(path, "wb") as handle:
            handle.write(tomli_w.dumps(payload).encode("utf-8"))

    @classmethod
    def from_toml(cls, path: str) -> "GenomicFigure":
        """Load a `GenomicFigure` definition from a TOML file."""
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

        track_registry = cls._track_registry()
        tracks_payload = payload.get("tracks", [])

        if isinstance(tracks_payload, list):
            for track_spec in tracks_payload:
                type_name = track_spec["type"]
                params = track_spec.get("params", {})
                if type_name not in track_registry:
                    raise ValueError(f"Unknown track type in TOML: {type_name}")
                fig.add_track(track_registry[type_name](**params))
            return fig

        if isinstance(tracks_payload, dict):
            for type_name, raw_params in tracks_payload.items():
                if type_name not in track_registry:
                    raise ValueError(f"Unknown track type in TOML: {type_name}")

                params_list = raw_params if isinstance(raw_params, list) else [raw_params]
                for params in params_list:
                    fig.add_track(track_registry[type_name](**(params or {})))
            return fig

        raise ValueError("Invalid TOML format: 'tracks' must be a list or table")

    def save(
        self,
        path: str | Path,
        region: str | GenomicRegion,
        dpi: int = 600,
        **kwargs: Any,
    ) -> None:
        """Render one region and write it to disk."""
        fig = self.plot(region, show=False)
        if fig is not None:
            fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
            logger.info(f"Saved figure to {path}")
            plt.close(fig)

    def __repr__(self) -> str:
        parts = [f"tracks={len(self.tracks)}"]
        if self.width != 12:
            parts.append(f"width={self.width}")
        if self.theme is None:
            parts.append("theme=none")
        elif self.theme == Theme.minimal():
            parts.append("theme=minimal")
        elif self.theme == Theme.publication():
            parts.append("theme=publication")
        elif self.theme == Theme.default():
            parts.append("theme=default")
        else:
            parts.append("theme=custom")
        parts.append(f"autoscale={self._autoscale}")
        parts.append(f"highlights={len(self._highlight_regions)}")
        return f"GenomicFigure({', '.join(parts)})"

    def _repr_html_(self) -> str:
        import html

        if not self.tracks:
            return (
                "<div><em>No tracks added yet. Use .bigwig(), .genes(), .axis(), etc.</em></div>"
            )

        rows = []
        for index, track in enumerate(self.tracks):
            title = getattr(track, "title", "") or ""
            data_value = getattr(track, "data", None)
            data_display = "" if data_value is None else str(data_value)
            rows.append(
                "<tr>"
                f"<td>{index}</td>"
                f"<td>{html.escape(track.__class__.__name__)}</td>"
                f"<td>{html.escape(str(title))}</td>"
                f"<td>{html.escape(data_display)}</td>"
                f"<td>{html.escape(str(getattr(track, 'height', '')))}</td>"
                "</tr>"
            )

        table = (
            "<table>"
            "<thead><tr><th>#</th><th>Type</th><th>Title</th><th>Data</th><th>Height</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody>"
            "</table>"
        )
        return table
