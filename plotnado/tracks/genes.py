"""
Genes track for displaying gene annotations.
"""

import importlib.resources
import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import matplotlib.axes
import matplotlib.font_manager
import matplotlib.markers
import matplotlib.patches
import matplotlib.textpath
import pandas as pd
from pydantic import ConfigDict, Field, BaseModel

from .base import Track
from .enums import DisplayMode, GeneLabelOverlapStrategy, GeneLabelStyle, PlotStyle, TrackType
from .region import GenomicRegion
from .utils import clean_axis, read_bed_regions, read_gtf_regions
from .aesthetics import BaseAesthetics
from .registry import registry


@dataclass
class LabelPlacement:
    text: str
    x: float
    y: float
    ha: str
    visible: bool
    row_index: int
    connector_x: float
    connector_y: float
    draw_connector: bool


class GenesAesthetics(BaseAesthetics):
    """Aesthetics configuration for gene tracks.

    Inherits color, alpha, and linewidth from BaseAesthetics.
    """

    style: PlotStyle = Field(default=PlotStyle.STD, description="Gene rendering style preset.")
    fill: bool = Field(default=True, description="Fill exon bodies instead of outlines only.")
    display: DisplayMode = Field(
        default=DisplayMode.COLLAPSED,
        description="Collapsed uses one row; expanded stacks overlapping genes.",
    )
    minimum_gene_length: int = Field(
        default=0,
        description="Minimum feature length (bp) required to draw a gene.",
    )
    max_number_of_rows: int = Field(
        default=4,
        description="Maximum rows shown when display mode is expanded.",
    )
    interval_height: float = Field(default=0.1, description="Vertical thickness of exon rectangles.")
    arrow_size: float = Field(default=0.12, description="Arrow marker size for transcript direction cues.")
    exon_linewidth: float = Field(default=0.8, description="Line width for exon outlines.")
    exon_edge_color: str = Field(default="black", description="Edge color for exon rectangles.")
    exon_color: str = Field(default="black", description="Fill color for exon rectangles.")
    intron_linewidth: float = Field(default=0.5, description="Line width for intron connector lines.")
    intron_color: str = Field(default="black", description="Color for intron connector lines.")
    chevron_height_ratio: float = Field(
        default=0.22,
        description="Chevron height relative to track row height.",
    )
    chevron_vertical_offset_ratio: float = Field(
        default=0.0,
        description="Vertical offset ratio used when drawing directional chevrons.",
    )
    chevron_width_fraction: float = Field(
        default=0.035,
        description="Default chevron width as fraction of viewport width.",
    )
    chevron_min_width_bp: float = Field(
        default=40.0,
        description="Minimum chevron width in base pairs.",
    )
    chevron_margin_bp: float = Field(
        default=20.0,
        description="Margin near transcript edges where chevrons are omitted.",
    )
    chevron_target_spacing_bp: float = Field(
        default=9000.0,
        description="Target spacing in bp between directional chevrons.",
    )
    chevron_max_count: int = Field(default=10, description="Maximum number of chevrons drawn per gene body.")
    gene_label_font_size: int = Field(default=10, description="Font size for gene name labels.")
    gene_label_style: GeneLabelStyle = Field(
        default=GeneLabelStyle.ITALIC,
        description="Font style for gene name labels.",
    )
    show_labels: bool = Field(default=True, description="Draw gene name labels.")
    label_overlap_strategy: GeneLabelOverlapStrategy = Field(
        default=GeneLabelOverlapStrategy.AUTO,
        description=(
            "How to handle overlapping labels: auto resolves by display mode, smart uses "
            "lane stacking + horizontal nudging, stagger alternates vertical offsets, "
            "suppress hides collisions, auto_expand switches collapsed mode to expanded."
        ),
    )
    label_max_chars: int = Field(
        default=20,
        description="Maximum characters to show for each gene label before truncation.",
    )
    label_offset_fraction: float = Field(
        default=0.005,
        description="Horizontal gene label offset as a fraction of viewport width.",
    )
    label_stagger_offset: float = Field(
        default=0.15,
        description="Vertical offset (in row units) used by staggered label placement.",
    )
    label_vertical_offset: float = Field(
        default=0.14,
        description="Vertical offset (in row units) for non-stagger label placement.",
    )
    label_min_overlap_bp: int = Field(
        default=200,
        description=(
            "Do not draw labels for genes clipped at viewport edges when visible overlap is below this bp threshold."
        ),
    )
    label_min_overlap_fraction: float = Field(
        default=0.15,
        description=(
            "Do not draw labels for clipped edge genes when overlap is below this fraction of total gene length."
        ),
    )
    label_connectors: bool = Field(
        default=False,
        description="Draw connector lines between displaced labels and their gene models.",
    )
    label_connector_linewidth: float = Field(
        default=0.7,
        description="Line width for label connector lines.",
    )

    model_config = ConfigDict(use_enum_values=True)


@registry.register(TrackType.GENE, aliases=["genes"])
class Genes(Track):
    """Track for displaying gene annotations from BED12 or GTF files."""

    data: Path | str | pd.DataFrame | None = Field(
        default=None,
        description="Gene annotation source (BED12/GTF file path or DataFrame).",
    )
    genome: str | None = Field(
        default=None,
        description="Bundled genome key used when data is not provided (e.g., hg38).",
    )
    aesthetics: GenesAesthetics = Field(
        default_factory=GenesAesthetics,
        description="Visual styling options for genes, exons, introns, and labels.",
    )
    height: float = Field(default=1.5, description="Relative panel height for this track.")

    row_scale: float = Field(default=1.0, description="Internal row scaling factor used during layout.")
    small_relative: float = Field(
        default=0.01,
        description="Internal small-offset ratio used for annotation geometry.",
    )
    gene_count: int = Field(default=0, description="Number of genes drawn in the last plotting pass.")

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @staticmethod
    def _enum_token(value: object) -> str:
        if isinstance(value, Enum):
            return str(value.value).lower()
        if hasattr(value, "value"):
            raw_value = getattr(value, "value")
            return str(raw_value).lower()
        raw = str(value)
        if "." in raw:
            raw = raw.split(".")[-1]
        return raw.lower()

    def _fetch_genes_from_package(self, gr: GenomicRegion) -> pd.DataFrame:
        try:
            bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
            mapping_path = bed_prefix / "genes.json"
            with open(mapping_path) as handle:
                gene_files = json.load(handle)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                "The genes database is not available. Please provide a data file."
            ) from exc

        if self.genome not in gene_files:
            raise ValueError(
                f"Genome {self.genome} not found in the genes database. "
                f"Available genomes: {list(gene_files.keys())}"
            )

        gene_file = bed_prefix / gene_files[self.genome]
        return self._fetch_from_disk_bed12(gr, Path(gene_file))

    @staticmethod
    def _parse_int_list(value: object) -> list[int]:
        if isinstance(value, list):
            return [int(v) for v in value]
        if pd.isna(value):
            return []
        raw = str(value).strip(",")
        if not raw:
            return []
        return [int(v) for v in raw.split(",") if str(v).strip()]

    def _fetch_from_disk_bed12(
        self, gr: GenomicRegion, file_path: Path | None = None
    ) -> pd.DataFrame:
        df = read_bed_regions(
            str(file_path or self.data), gr.chromosome, gr.start, gr.end
        )

        if df.empty:
            return pd.DataFrame(
                columns=[
                    "chrom",
                    "start",
                    "end",
                    "geneid",
                    "strand",
                    "block_count",
                    "block_sizes",
                    "block_starts",
                ]
            )

        # Normalize BED12 naming variants from different parsers (pyranges, BigBed).
        rename_map = {
            "name": "geneid",
            "Name": "geneid",
            "field_1": "geneid",
            "score": "score",
            "Score": "score",
            "field_2": "score",
            "Strand": "strand",
            "field_3": "strand",
            "itemRgb": "rgb",
            "ItemRGB": "rgb",
            "field_6": "rgb",
            "blockCount": "block_count",
            "BlockCount": "block_count",
            "field_7": "block_count",
            "blockSizes": "block_sizes",
            "BlockSizes": "block_sizes",
            "field_8": "block_sizes",
            "blockStarts": "block_starts",
            "BlockStarts": "block_starts",
            "field_9": "block_starts",
            "thickStart": "thick_start",
            "ThickStart": "thick_start",
            "field_4": "thick_start",
            "thickEnd": "thick_end",
            "ThickEnd": "thick_end",
            "field_5": "thick_end",
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        if "geneid" not in df.columns:
            df["geneid"] = "gene"
        if "strand" not in df.columns:
            df["strand"] = "+"
        if "block_sizes" not in df.columns:
            df["block_sizes"] = (df["end"] - df["start"]).astype(int).astype(str)
        if "block_starts" not in df.columns:
            df["block_starts"] = "0"

        df["block_starts"] = df["block_starts"].apply(self._parse_int_list)
        df["block_sizes"] = df["block_sizes"].apply(self._parse_int_list)

        if "block_count" in df.columns:
            df["block_count"] = (
                pd.to_numeric(df["block_count"], errors="coerce")
                .fillna(df["block_starts"].apply(lambda values: len(values) or 1))
                .astype(int)
            )
        else:
            df["block_count"] = df["block_starts"].apply(lambda values: len(values) or 1)
        return df

    def _fetch_from_disk_gtf(self, gr: GenomicRegion) -> pd.DataFrame:
        gtf_df = read_gtf_regions(str(self.data), gr.chromosome, gr.start, gr.end)
        if gtf_df.empty:
            return pd.DataFrame(
                columns=[
                    "chrom",
                    "start",
                    "end",
                    "geneid",
                    "strand",
                    "block_count",
                    "block_sizes",
                    "block_starts",
                ]
            )

        if "Feature" in gtf_df.columns:
            gtf_df = gtf_df.loc[gtf_df["Feature"].astype(str).str.lower() == "exon"].copy()
        elif "feature" in gtf_df.columns:
            gtf_df = gtf_df.loc[gtf_df["feature"].astype(str).str.lower() == "exon"].copy()

        if gtf_df.empty:
            return pd.DataFrame(
                columns=[
                    "chrom",
                    "start",
                    "end",
                    "geneid",
                    "strand",
                    "block_count",
                    "block_sizes",
                    "block_starts",
                ]
            )

        gene_column = "gene_id" if "gene_id" in gtf_df.columns else "geneid"
        strand_column = "Strand" if "Strand" in gtf_df.columns else "strand"
        start_column = "start" if "start" in gtf_df.columns else "Start"
        end_column = "end" if "end" in gtf_df.columns else "End"
        chrom_column = "chrom" if "chrom" in gtf_df.columns else "Chromosome"

        records: list[dict] = []
        for gene_id, group in gtf_df.groupby(gene_column):
            group = group.sort_values(start_column)
            gene_start = int(group[start_column].min())
            gene_end = int(group[end_column].max())
            block_starts = (group[start_column].astype(int) - gene_start).tolist()
            block_sizes = (
                group[end_column].astype(int) - group[start_column].astype(int)
            ).tolist()
            records.append(
                {
                    "chrom": group[chrom_column].iloc[0],
                    "start": gene_start,
                    "end": gene_end,
                    "geneid": str(gene_id),
                    "strand": group[strand_column].iloc[0] if strand_column in group.columns else "+",
                    "block_count": len(block_starts),
                    "block_starts": block_starts,
                    "block_sizes": block_sizes,
                }
            )

        return pd.DataFrame(records)

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:
        if self.data is None and self.genome is None:
            raise ValueError("Either data or genome must be provided")

        if self.data is None:
            data = self._fetch_genes_from_package(gr)
        elif isinstance(self.data, pd.DataFrame):
            data = self.data.copy()
        else:
            data_str = str(self.data).lower()
            if data_str.endswith(".gtf") or data_str.endswith(".gtf.gz"):
                data = self._fetch_from_disk_gtf(gr)
            elif (
                data_str.endswith(".bed")
                or data_str.endswith(".bed.gz")
                or data_str.endswith(".bb")
                or data_str.endswith(".bigbed")
            ):
                data = self._fetch_from_disk_bed12(gr)
            else:
                raise ValueError(
                    "Unsupported file format. Only BED/BED.GZ/BigBed and GTF files are supported."
                )

        if self.minimum_gene_length > 0 and not data.empty:
            data = data.query(f"end - start >= {self.minimum_gene_length}")
        return data

    def _allocate_row_index(
        self, row_last_positions: list[int], start_bp: int, end_bp: int
    ) -> int:
        for index, last_end in enumerate(row_last_positions):
            if last_end < start_bp:
                row_last_positions[index] = end_bp
                return index
        row_last_positions.append(end_bp)
        return len(row_last_positions) - 1

    @staticmethod
    def _truncate_label(text: str, max_chars: int) -> str:
        if max_chars <= 0:
            return ""
        if len(text) <= max_chars:
            return text
        if max_chars == 1:
            return "…"
        return f"{text[: max_chars - 1]}…"

    def _estimate_bp_per_px(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> float:
        viewport_bp = max(float(gr.length), 1.0)
        axes_width_px = 1000.0
        try:
            figure = ax.figure
            fig_size_inches = figure.get_size_inches()
            fig_width_px = max(float(fig_size_inches[0] * figure.dpi), 1.0)
            axes_fraction = max(float(ax.get_position().width), 1e-4)
            axes_width_px = max(fig_width_px * axes_fraction, 1.0)
            if axes_width_px < 10.0:
                axes_width_px = 1000.0
        except (AttributeError, TypeError, ValueError, IndexError):
            axes_width_px = 1000.0
        return viewport_bp / axes_width_px

    @staticmethod
    def _resolve_figure_dpi(ax: matplotlib.axes.Axes) -> float:
        try:
            return max(float(ax.figure.dpi), 72.0)
        except (AttributeError, TypeError, ValueError):
            return 100.0

    def _measure_label_width_bp(
        self,
        text: str,
        bp_per_px: float,
        dpi: float,
    ) -> float:
        if not text:
            return 0.0

        fallback_px = max(float(self.gene_label_font_size) * 1.0 * len(text), 1.0)
        fallback_bp = fallback_px * bp_per_px

        try:
            font_prop = matplotlib.font_manager.FontProperties(
                size=float(self.gene_label_font_size),
                style=str(self.gene_label_style),
            )
            width_points = matplotlib.textpath.TextPath((0, 0), text, prop=font_prop).get_extents().width
            width_px = max(float(width_points) * (dpi / 72.0), 1.0)
            return max(width_px * bp_per_px, fallback_bp)
        except (TypeError, ValueError, RuntimeError):
            return fallback_bp

    def _resolve_effective_label_strategy(self, display_mode: str) -> str:
        strategy = self._enum_token(self.label_overlap_strategy)
        if strategy != "auto":
            return strategy
        if display_mode == "collapsed":
            return "smart"
        return "smart"

    @staticmethod
    def _label_interval_from_anchor(anchor_x: float, width: float, ha: str) -> tuple[float, float]:
        if ha == "left":
            return anchor_x, anchor_x + width
        if ha == "right":
            return anchor_x - width, anchor_x
        half_width = width / 2.0
        return anchor_x - half_width, anchor_x + half_width

    def _resolve_label_geometry(
        self,
        raw_label: str,
        anchor_x: float,
        ha: str,
        gr: GenomicRegion,
        bp_per_px: float,
        dpi: float,
    ) -> tuple[str, float, str, float, float, bool]:
        max_chars = max(1, int(self.label_max_chars))
        fallback_start, fallback_end = self._label_interval_from_anchor(anchor_x, 1.0, ha)

        for char_count in range(max_chars, 0, -1):
            label = self._truncate_label(raw_label, char_count)
            width = max(1.0, self._measure_label_width_bp(label, bp_per_px=bp_per_px, dpi=dpi))
            label_start, label_end = self._label_interval_from_anchor(anchor_x, width, ha)
            fallback_start, fallback_end = label_start, label_end
            if label_start >= gr.start and label_end <= gr.end:
                return label, anchor_x, ha, label_start, label_end, True

        return "", anchor_x, ha, fallback_start, fallback_end, False

    def _label_candidates(
        self,
        raw_label: str,
        gene_start: int,
        gene_end: int,
        gr: GenomicRegion,
        bp_per_px: float,
        dpi: float,
    ) -> list[tuple[str, float, str, float, float, bool]]:
        center_x = (float(gene_start) + float(gene_end)) / 2.0
        offset_bp = max(0.0, float(self.label_offset_fraction) * float(gr.length))
        right_anchor = min(float(gr.end), float(gene_end) + offset_bp)
        left_anchor = max(float(gr.start), float(gene_start) - offset_bp)

        candidates = [
            self._resolve_label_geometry(
                raw_label=raw_label,
                anchor_x=center_x,
                ha="center",
                gr=gr,
                bp_per_px=bp_per_px,
                dpi=dpi,
            ),
            self._resolve_label_geometry(
                raw_label=raw_label,
                anchor_x=right_anchor,
                ha="left",
                gr=gr,
                bp_per_px=bp_per_px,
                dpi=dpi,
            ),
            self._resolve_label_geometry(
                raw_label=raw_label,
                anchor_x=left_anchor,
                ha="right",
                gr=gr,
                bp_per_px=bp_per_px,
                dpi=dpi,
            ),
        ]
        return [candidate for candidate in candidates if candidate[5]]

    def _should_label_gene(self, gene: pd.Series, gr: GenomicRegion, bp_per_px: float) -> bool:
        gene_start = int(gene["start"])
        gene_end = int(gene["end"])
        gene_len = max(1, gene_end - gene_start)
        overlap_start = max(gene_start, int(gr.start))
        overlap_end = min(gene_end, int(gr.end))
        overlap_bp = max(0, overlap_end - overlap_start)
        overlap_fraction = overlap_bp / float(gene_len)
        clipped_at_edge = gene_start < int(gr.start) or gene_end > int(gr.end)

        # Skip labels for genes with no visible body in the current viewport.
        if overlap_bp < max(1.0, bp_per_px):
            return False

        if (
            clipped_at_edge
            and overlap_bp < int(self.label_min_overlap_bp)
            and overlap_fraction < float(self.label_min_overlap_fraction)
        ):
            return False
        return True

    def _compute_label_placements(
        self,
        genes_df: pd.DataFrame,
        row_indices: list[int],
        ax: matplotlib.axes.Axes,
        gr: GenomicRegion,
    ) -> tuple[list[LabelPlacement], bool]:
        placements: list[LabelPlacement] = []
        if not self.show_labels:
            return [
                LabelPlacement(
                    text="",
                    x=0.0,
                    y=0.0,
                    ha="left",
                    visible=False,
                    row_index=row,
                    connector_x=0.0,
                    connector_y=0.0,
                    draw_connector=False,
                )
                for row in row_indices
            ], False

        strategy = self._resolve_effective_label_strategy(self._enum_token(self.display))
        bp_per_px = self._estimate_bp_per_px(ax, gr)
        dpi = self._resolve_figure_dpi(ax)
        row_last_end: dict[int, float] = {}
        row_lane_last_end: dict[tuple[int, int], float] = {}
        row_next_lane: dict[int, int] = {}
        smart_lane_last_end: dict[tuple[int, int], float] = {}
        had_collisions = False
        min_gap_bp = max(2.0 * bp_per_px, 0.4 * float(self.gene_label_font_size) * bp_per_px)
        lane_step = max(0.01, float(self.label_stagger_offset) * float(self.row_scale))

        for (_, gene), row_index in zip(genes_df.iterrows(), row_indices):
            raw_label = str(gene.get("geneid", "gene"))
            base_y = self.get_y_pos(row_index)
            center_x = (float(gene["start"]) + float(gene["end"])) / 2.0
            if not self._should_label_gene(gene, gr, bp_per_px):
                placements.append(
                    LabelPlacement(
                        text="",
                        x=center_x,
                        y=base_y,
                        ha="center",
                        visible=False,
                        row_index=row_index,
                        connector_x=center_x,
                        connector_y=base_y,
                        draw_connector=False,
                    )
                )
                continue
            candidates = self._label_candidates(
                raw_label=raw_label,
                gene_start=int(gene["start"]),
                gene_end=int(gene["end"]),
                gr=gr,
                bp_per_px=bp_per_px,
                dpi=dpi,
            )

            if not candidates:
                placements.append(
                    LabelPlacement(
                        text="",
                        x=center_x,
                        y=base_y,
                        ha="center",
                        visible=False,
                        row_index=row_index,
                        connector_x=center_x,
                        connector_y=base_y,
                        draw_connector=False,
                    )
                )
                continue

            if strategy == "smart":
                chosen: LabelPlacement | None = None
                for lane_index in range(max(1, len(genes_df) + 1)):
                    y_pos = (
                        base_y
                        + (float(self.label_vertical_offset) * float(self.row_scale))
                        + (lane_index * lane_step)
                    )
                    lane_key = (row_index, lane_index)
                    last_end = smart_lane_last_end.get(lane_key, float("-inf"))

                    for text, x_pos, ha, label_start, label_end, _ in candidates:
                        if label_start <= (last_end + min_gap_bp):
                            continue

                        smart_lane_last_end[lane_key] = label_end
                        displaced = lane_index > 0
                        if displaced:
                            had_collisions = True
                        chosen = LabelPlacement(
                            text=text,
                            x=x_pos,
                            y=y_pos,
                            ha=ha,
                            visible=True,
                            row_index=row_index,
                            connector_x=center_x,
                            connector_y=base_y,
                            draw_connector=(
                                bool(self.label_connectors)
                                and (
                                    abs(x_pos - center_x) > 1e-6
                                    or abs(y_pos - base_y) > 1e-6
                                )
                            ),
                        )
                        break

                    if chosen is not None:
                        break

                if chosen is not None:
                    placements.append(chosen)
                    continue

                had_collisions = True
                placements.append(
                    LabelPlacement(
                        text=candidates[0][0],
                        x=candidates[0][1],
                        y=base_y,
                        ha=candidates[0][2],
                        visible=False,
                        row_index=row_index,
                        connector_x=center_x,
                        connector_y=base_y,
                        draw_connector=False,
                    )
                )
                continue

            text, x_pos, ha, label_start, label_end, _ = candidates[0]
            if strategy == "stagger":
                preferred_lane = row_next_lane.get(row_index, 0)
                row_next_lane[row_index] = 1 - preferred_lane
                lanes = [preferred_lane, 1 - preferred_lane]
                lane_used: int | None = None
                for lane in lanes:
                    lane_last_end = row_lane_last_end.get((row_index, lane), float("-inf"))
                    if label_start > lane_last_end:
                        lane_used = lane
                        row_lane_last_end[(row_index, lane)] = label_end
                        break

                if lane_used is None:
                    had_collisions = True
                    placements.append(
                        LabelPlacement(
                            text=text,
                            x=x_pos,
                            y=base_y,
                            ha=ha,
                            visible=False,
                            row_index=row_index,
                            connector_x=center_x,
                            connector_y=base_y,
                            draw_connector=False,
                        )
                    )
                    continue

                y_offset = float(self.label_stagger_offset) * float(self.row_scale)
                y_pos = base_y + y_offset if lane_used == 0 else base_y - y_offset
                placements.append(
                    LabelPlacement(
                        text=text,
                        x=x_pos,
                        y=y_pos,
                        ha=ha,
                        visible=True,
                        row_index=row_index,
                        connector_x=center_x,
                        connector_y=base_y,
                        draw_connector=(
                            bool(self.label_connectors)
                            and (
                                abs(x_pos - center_x) > 1e-6
                                or abs(y_pos - base_y) > 1e-6
                            )
                        ),
                    )
                )
                continue

            last_end = row_last_end.get(row_index, float("-inf"))
            overlaps = label_start <= last_end
            if overlaps:
                had_collisions = True
                if strategy == "suppress" or strategy == "auto_expand":
                    placements.append(
                        LabelPlacement(
                            text=text,
                            x=x_pos,
                            y=base_y,
                            ha=ha,
                            visible=False,
                            row_index=row_index,
                            connector_x=center_x,
                            connector_y=base_y,
                            draw_connector=False,
                        )
                    )
                    continue

            row_last_end[row_index] = label_end
            y_offset = float(self.label_vertical_offset) * float(self.row_scale)
            y_pos = base_y + y_offset
            placements.append(
                LabelPlacement(
                    text=text,
                    x=x_pos,
                    y=y_pos,
                    ha=ha,
                    visible=True,
                    row_index=row_index,
                    connector_x=center_x,
                    connector_y=base_y,
                    draw_connector=(
                        bool(self.label_connectors)
                        and (
                            abs(x_pos - center_x) > 1e-6
                            or abs(y_pos - base_y) > 1e-6
                        )
                    ),
                )
            )

        return placements, had_collisions

    def _label_footprint_end(
        self,
        gene: pd.Series,
        gr: GenomicRegion,
        bp_per_px: float,
        dpi: float,
    ) -> int:
        if not self.show_labels:
            return int(gene["end"])
        raw_label = str(gene.get("geneid", "gene"))
        label = self._truncate_label(raw_label, max(1, int(self.label_max_chars)))
        label_width_bp = max(1.0, self._measure_label_width_bp(label, bp_per_px=bp_per_px, dpi=dpi))
        center_x = (float(gene["start"]) + float(gene["end"])) / 2.0
        label_end = center_x + (label_width_bp / 2.0)
        return int(max(float(gene["end"]), label_end))

    def get_y_pos(self, row_index: int) -> float:
        return row_index * self.row_scale + self.row_scale / 2

    def _draw_gene_with_introns(
        self, ax: matplotlib.axes.Axes, gene: pd.Series, ypos: float
    ) -> None:
        block_starts = gene.block_starts if isinstance(gene.block_starts, list) else []
        block_sizes = gene.block_sizes if isinstance(gene.block_sizes, list) else []
        if not block_starts or not block_sizes:
            self._draw_gene_simple(ax, gene, ypos)
            return

        for index, (start, size) in enumerate(zip(block_starts, block_sizes)):
            block_start = gene.start + start
            rect = matplotlib.patches.Rectangle(
                (block_start, ypos - self.interval_height / 2),
                size,
                self.interval_height,
                linewidth=self.exon_linewidth,
                edgecolor=self.exon_edge_color,
                facecolor=self.exon_color,
                alpha=self.alpha,
                zorder=2,
            )
            ax.add_patch(rect)

            if index < len(block_starts) - 1:
                intron_start = block_start + size
                intron_end = gene.start + block_starts[index + 1]
                ax.plot(
                    [intron_start, intron_end],
                    [ypos, ypos],
                    color=self.intron_color,
                    linewidth=self.intron_linewidth,
                    zorder=1,
                )
                self._draw_intron_chevrons(
                    ax=ax,
                    intron_start=intron_start,
                    intron_end=intron_end,
                    ypos=ypos,
                    strand=str(gene.get("strand", "+")),
                )

    def _draw_intron_chevrons(
        self,
        ax: matplotlib.axes.Axes,
        intron_start: float,
        intron_end: float,
        ypos: float,
        strand: str,
    ) -> None:
        intron_len = intron_end - intron_start
        if intron_len <= 0:
            return

        chevron_height = self._resolve_chevron_height(ax=ax)
        chevron_y = ypos - (chevron_height * self.chevron_vertical_offset_ratio)
        direction = self._strand_direction(strand)
        marker = self._strand_marker(direction)
        chevron_width = self._resolve_chevron_width(
            ax=ax,
            intron_len=intron_len,
            chevron_height=chevron_height,
        )
        marker_size = self._resolve_chevron_marker_size(ax=ax, chevron_height=chevron_height)
        marker_edge_width = self._resolve_chevron_marker_edge_width()
        marker_left_extent, marker_right_extent = self._resolve_chevron_marker_extents(
            ax=ax,
            marker=marker,
            marker_size=marker_size,
            marker_edge_width=marker_edge_width,
        )
        boundary_padding = self._resolve_chevron_boundary_padding(
            ax=ax,
            marker_edge_width=marker_edge_width,
        )
        left_extent = max(chevron_width / 2, marker_left_extent)
        right_extent = max(chevron_width / 2, marker_right_extent)
        allowed_left = intron_start + boundary_padding
        allowed_right = intron_end - boundary_padding
        start_center = allowed_left + self.chevron_margin_bp + left_extent
        end_center = allowed_right - self.chevron_margin_bp - right_extent
        usable_len = end_center - start_center
        if usable_len <= 0:
            return

        target_spacing_bp = max(1.0, self.chevron_target_spacing_bp)
        max_count = max(1, self.chevron_max_count)
        chevron_count = max(1, min(max_count, int(usable_len / target_spacing_bp) + 1))
        if chevron_count > 1 and chevron_count % 2 == 0:
            if chevron_count < max_count:
                chevron_count += 1
            else:
                chevron_count -= 1

        midpoint = (intron_start + intron_end) / 2.0

        if chevron_count == 1:
            centers = [midpoint]
        else:
            max_symmetric_spacing = (end_center - start_center) / (chevron_count - 1)
            spacing = max_symmetric_spacing
            half_span = spacing * ((chevron_count - 1) / 2.0)
            first_center = midpoint - half_span
            centers = [first_center + (index * spacing) for index in range(chevron_count)]

        clip_rect = self._build_chevron_clip_rect(
            ax=ax,
            intron_start=allowed_left,
            intron_end=allowed_right,
            chevron_y=chevron_y,
            chevron_height=chevron_height,
        )
        for center in centers:
            if (
                center - marker_left_extent < allowed_left
                or center + marker_right_extent > allowed_right
            ):
                continue
            line_artists = ax.plot(
                [center],
                [chevron_y],
                linestyle="None",
                marker=marker,
                markersize=marker_size,
                markerfacecolor="none",
                markeredgecolor=self.intron_color,
                markeredgewidth=marker_edge_width,
                zorder=1.5,
            )
            if isinstance(line_artists, (list, tuple)) and line_artists:
                line_artists[0].set_clip_path(clip_rect.get_path(), clip_rect.get_transform())

    def _resolve_chevron_width(
        self,
        ax: matplotlib.axes.Axes,
        intron_len: float,
        chevron_height: float,
    ) -> float:
        width_from_span = max(self.chevron_min_width_bp, intron_len * self.chevron_width_fraction)
        width_from_display = self._chevron_width_from_display(ax=ax, chevron_height=chevron_height)
        if width_from_display is None:
            return width_from_span
        return min(width_from_span, width_from_display)

    def _resolve_chevron_marker_size(
        self,
        ax: matplotlib.axes.Axes,
        chevron_height: float,
    ) -> float:
        scales = self._data_pixel_scales(ax)
        dpi = float(getattr(getattr(ax, "figure", None), "dpi", 72.0) or 72.0)
        if scales is None:
            return max(3.5, chevron_height * 54.0)

        _, px_per_y = scales
        marker_height_px = max(6.0, 1.6 * chevron_height * px_per_y)
        return marker_height_px * 72.0 / dpi

    def _resolve_chevron_marker_edge_width(self) -> float:
        return max(0.4, 0.8 * self.intron_linewidth)

    def _resolve_chevron_marker_extents(
        self,
        ax: matplotlib.axes.Axes,
        marker: int,
        marker_size: float,
        marker_edge_width: float,
    ) -> tuple[float, float]:
        scales = self._data_pixel_scales(ax)
        if scales is None:
            fallback = self.chevron_min_width_bp / 2
            return fallback, fallback

        px_per_x, _ = scales
        dpi = float(getattr(getattr(ax, "figure", None), "dpi", 72.0) or 72.0)
        marker_style = matplotlib.markers.MarkerStyle(marker)
        marker_path = marker_style.get_path().transformed(marker_style.get_transform())
        vertices = marker_path.vertices
        min_x = abs(float(vertices[:, 0].min()))
        max_x = abs(float(vertices[:, 0].max()))
        points_to_px = dpi / 72.0
        left_extent_px = (min_x * marker_size * points_to_px) + (marker_edge_width / 2.0)
        right_extent_px = (max_x * marker_size * points_to_px) + (marker_edge_width / 2.0)
        return left_extent_px / px_per_x, right_extent_px / px_per_x

    def _resolve_chevron_boundary_padding(
        self,
        ax: matplotlib.axes.Axes,
        marker_edge_width: float,
    ) -> float:
        scales = self._data_pixel_scales(ax)
        if scales is None:
            return 0.0

        px_per_x, _ = scales
        padding_px = max(1.5, marker_edge_width)
        return padding_px / px_per_x

    def _build_chevron_clip_rect(
        self,
        ax: matplotlib.axes.Axes,
        intron_start: float,
        intron_end: float,
        chevron_y: float,
        chevron_height: float,
    ) -> matplotlib.patches.Rectangle:
        clip_half_height = max(self.row_scale, chevron_height * 4.0)
        return matplotlib.patches.Rectangle(
            (intron_start, chevron_y - clip_half_height),
            intron_end - intron_start,
            clip_half_height * 2.0,
            transform=ax.transData,
        )

    def _resolve_chevron_height(self, ax: matplotlib.axes.Axes) -> float:
        height_from_row = max(self.interval_height * self.chevron_height_ratio, 0.01)
        height_from_display = self._chevron_height_from_display(ax=ax)
        if height_from_display is None:
            return height_from_row
        return max(height_from_row, height_from_display)

    def _chevron_width_from_display(
        self,
        ax: matplotlib.axes.Axes,
        chevron_height: float,
    ) -> float | None:
        scales = self._data_pixel_scales(ax)
        if scales is None:
            return None
        px_per_x, px_per_y = scales

        chevron_height_px = chevron_height * px_per_y
        # Keep reduced-height chevrons visibly open instead of collapsing into a single stroke.
        min_visible_width_px = 8.0
        full_width_px = max(min_visible_width_px, 3.0 * chevron_height_px)
        return full_width_px / px_per_x

    def _chevron_height_from_display(
        self,
        ax: matplotlib.axes.Axes,
    ) -> float | None:
        scales = self._data_pixel_scales(ax)
        if scales is None:
            return None
        _, px_per_y = scales
        min_visible_half_height_px = max(3.0, 2.0 * self.intron_linewidth)
        return min_visible_half_height_px / px_per_y

    def _data_pixel_scales(self, ax: matplotlib.axes.Axes) -> tuple[float, float] | None:
        trans_data = getattr(ax, "transData", None)
        if trans_data is None:
            return None

        try:
            origin = trans_data.transform((0.0, 0.0))
            x_unit = trans_data.transform((1.0, 0.0))
            y_unit = trans_data.transform((0.0, 1.0))
        except Exception:
            return None

        px_per_x = abs(float(x_unit[0]) - float(origin[0]))
        px_per_y = abs(float(y_unit[1]) - float(origin[1]))
        if px_per_x <= 0 or px_per_y <= 0:
            return None
        return px_per_x, px_per_y

    @staticmethod
    def _strand_marker(direction: int) -> int:
        return (
            matplotlib.markers.CARETRIGHTBASE
            if direction >= 0
            else matplotlib.markers.CARETLEFTBASE
        )

    @staticmethod
    def _strand_direction(strand: str) -> int:
        strand_value = str(strand).strip().lower()
        if strand_value in {"-", "-1", "minus", "negative", "reverse", "rev"}:
            return -1
        return 1

    def _draw_gene_simple(
        self, ax: matplotlib.axes.Axes, gene: pd.Series, ypos: float
    ) -> None:
        rect = matplotlib.patches.Rectangle(
            (gene.start, ypos - self.interval_height / 2),
            gene.end - gene.start,
            self.interval_height,
            linewidth=self.exon_linewidth,
            edgecolor=self.exon_edge_color,
            facecolor=self.exon_color,
            alpha=self.alpha,
        )
        ax.add_patch(rect)

    def _draw_gene_feature(
        self,
        ax: matplotlib.axes.Axes,
        gene: pd.Series,
        row_index: int,
        label: LabelPlacement | None = None,
    ) -> None:
        ypos = self.get_y_pos(row_index)
        block_count = int(gene.get("block_count", 1))
        if block_count > 1:
            self._draw_gene_with_introns(ax, gene, ypos)
        else:
            self._draw_gene_simple(ax, gene, ypos)

        if not self.show_labels:
            return

        if label is not None:
            if not label.visible:
                return
            if label.draw_connector and self.label_connectors:
                ax.plot(
                    [label.connector_x, label.x],
                    [label.connector_y, label.y],
                    color=self.intron_color,
                    linewidth=self.label_connector_linewidth,
                    zorder=0.5,
                )
            ax.text(
                label.x,
                label.y,
                label.text,
                va="center",
                ha=label.ha,
                fontsize=self.gene_label_font_size,
                style=self.gene_label_style,
                color="black",
                clip_on=False,
            )
            return

        center_x = (float(gene["start"]) + float(gene["end"])) / 2.0
        ax.text(
            center_x,
            ypos + (float(self.label_vertical_offset) * float(self.row_scale)),
            str(gene.get("geneid", "gene")),
            va="center",
            ha="center",
            fontsize=self.gene_label_font_size,
            style=self.gene_label_style,
            color="black",
            clip_on=False,
        )

    def plot_genes(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        genes_df = self.fetch_data(gr)

        if genes_df.empty:
            ax.set_ylim(0, 1)
            ax.set_xlim(gr.start, gr.end)
            return

        genes_df = genes_df.sort_values(["start", "end"]).reset_index(drop=True)

        display_mode = self._enum_token(self.display)
        strategy = self._resolve_effective_label_strategy(display_mode)
        bp_per_px = self._estimate_bp_per_px(ax, gr)
        dpi = self._resolve_figure_dpi(ax)

        def layout_rows(mode: str) -> tuple[list[int], int]:
            row_last_positions: list[int] = []
            row_indices: list[int] = []
            max_row_index = 0
            for _, gene in genes_df.iterrows():
                if mode == "expanded":
                    end_bp = int(gene["end"])
                    if self.show_labels:
                        end_bp = self._label_footprint_end(gene, gr, bp_per_px, dpi)
                    row_index = self._allocate_row_index(
                        row_last_positions,
                        int(gene["start"]),
                        end_bp,
                    )
                    if strategy != "smart":
                        row_index = min(row_index, self.max_number_of_rows - 1)
                else:
                    row_index = 0
                row_indices.append(row_index)
                max_row_index = max(max_row_index, row_index)
            return row_indices, max_row_index

        row_indices, max_row_index = layout_rows(display_mode)
        placements, had_collisions = self._compute_label_placements(genes_df, row_indices, ax, gr)

        if (
            strategy == "auto_expand"
            and display_mode == "collapsed"
            and had_collisions
        ):
            row_indices, max_row_index = layout_rows("expanded")
            placements, _ = self._compute_label_placements(genes_df, row_indices, ax, gr)

        for ((_, gene), row_index, placement) in zip(
            genes_df.iterrows(), row_indices, placements
        ):
            self._draw_gene_feature(ax, gene, row_index, placement)

        ax.set_xlim(gr.start, gr.end)
        rows = max(1, max_row_index + 1)
        label_top = rows * self.row_scale
        for placement in placements:
            if placement.visible:
                label_top = max(label_top, placement.y)
        ax.set_ylim(-0.1, label_top + 0.15)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        self.plot_genes(ax, gr)
        clean_axis(ax)
