"""
Genes track for displaying gene annotations.
"""

import importlib.resources
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.axes
import matplotlib.font_manager
import matplotlib.patches
import matplotlib.textpath
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from .base import Track
from .enums import DisplayMode, GeneLabelOverlapStrategy, GeneLabelStyle, PlotStyle
from .region import GenomicRegion
from .utils import clean_axis, read_bed_regions, read_gtf_regions


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


class GenesAesthetics(BaseModel):
    """Aesthetics configuration for gene tracks."""

    style: PlotStyle = Field(default=PlotStyle.STD, description="Gene rendering style preset.")
    color: str = Field(default="black", description="Primary color used for gene glyphs.")
    fill: bool = Field(default=True, description="Fill exon bodies instead of outlines only.")
    alpha: float = Field(default=1.0, description="Opacity for rendered gene bodies (0-1).")
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
    intron_linewidth: float = Field(default=0.8, description="Line width for intron connector lines.")
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
        default=6000.0,
        description="Target spacing in bp between directional chevrons.",
    )
    chevron_max_count: int = Field(default=14, description="Maximum number of chevrons drawn per gene body.")
    gene_label_font_size: int = Field(default=8, description="Font size for gene name labels.")
    gene_label_style: GeneLabelStyle = Field(
        default=GeneLabelStyle.ITALIC,
        description="Font style for gene name labels.",
    )
    show_labels: bool = Field(default=True, description="Draw gene name labels.")
    label_overlap_strategy: GeneLabelOverlapStrategy = Field(
        default=GeneLabelOverlapStrategy.STAGGER,
        description=(
            "How to handle overlapping labels: stagger alternates vertical offsets, "
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
    label_connectors: bool = Field(
        default=False,
        description="Draw connector lines between displaced labels and their gene models.",
    )
    label_connector_linewidth: float = Field(
        default=0.7,
        description="Line width for label connector lines.",
    )

    model_config = ConfigDict(use_enum_values=True)


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

    def _resolve_label_geometry(
        self,
        raw_label: str,
        gene_start: int,
        gene_end: int,
        gr: GenomicRegion,
        bp_per_px: float,
        dpi: float,
    ) -> tuple[str, float, str, float, float, bool]:
        max_chars = max(1, int(self.label_max_chars))

        def label_width_bp(text: str) -> float:
            return max(1.0, self._measure_label_width_bp(text, bp_per_px=bp_per_px, dpi=dpi))

        label = self._truncate_label(raw_label, max_chars)
        center_x = (float(gene_start) + float(gene_end)) / 2.0
        width = label_width_bp(label)
        label_start = center_x - (width / 2.0)
        label_end = center_x + (width / 2.0)
        if label_start >= gr.start and label_end <= gr.end:
            return label, center_x, "center", label_start, label_end, True

        available_width = max(0.0, min(2.0 * (center_x - gr.start), 2.0 * (gr.end - center_x)))
        if available_width <= 0:
            return "", center_x, "center", label_start, label_end, False

        avg_char_bp = max(self._measure_label_width_bp("M", bp_per_px=bp_per_px, dpi=dpi), 1.0)
        fit_chars = max(2, min(max_chars, int(available_width / max(avg_char_bp, 1e-6))))

        fit_label = self._truncate_label(raw_label, fit_chars)
        fit_width = label_width_bp(fit_label)
        while fit_chars > 1 and fit_width > available_width:
            fit_chars -= 1
            fit_label = self._truncate_label(raw_label, fit_chars)
            fit_width = label_width_bp(fit_label)

        if fit_chars <= 1 or fit_width > available_width:
            return "", center_x, "center", label_start, label_end, False

        fit_start = center_x - (fit_width / 2.0)
        fit_end = center_x + (fit_width / 2.0)
        return fit_label, center_x, "center", fit_start, fit_end, True

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

        strategy = str(self.label_overlap_strategy)
        bp_per_px = self._estimate_bp_per_px(ax, gr)
        dpi = self._resolve_figure_dpi(ax)
        row_last_end: dict[int, float] = {}
        row_lane_last_end: dict[tuple[int, int], float] = {}
        row_next_lane: dict[int, int] = {}
        had_collisions = False

        for (_, gene), row_index in zip(genes_df.iterrows(), row_indices):
            raw_label = str(gene.get("geneid", "gene"))
            base_y = self.get_y_pos(row_index)
            text, x_pos, ha, label_start, label_end, visible = self._resolve_label_geometry(
                raw_label=raw_label,
                gene_start=int(gene["start"]),
                gene_end=int(gene["end"]),
                gr=gr,
                bp_per_px=bp_per_px,
                dpi=dpi,
            )

            if not visible:
                placements.append(
                    LabelPlacement(
                        text="",
                        x=x_pos,
                        y=base_y,
                        ha=ha,
                        visible=False,
                        row_index=row_index,
                        connector_x=(float(gene["start"]) + float(gene["end"])) / 2.0,
                        connector_y=base_y,
                        draw_connector=False,
                    )
                )
                continue

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
                            connector_x=(float(gene["start"]) + float(gene["end"])) / 2.0,
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
                        connector_x=(float(gene["start"]) + float(gene["end"])) / 2.0,
                        connector_y=base_y,
                        draw_connector=False,
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
                            connector_x=(float(gene["start"]) + float(gene["end"])) / 2.0,
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
                    connector_x=(float(gene["start"]) + float(gene["end"])) / 2.0,
                    connector_y=base_y,
                    draw_connector=False,
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

        chevron_height = max(self.interval_height * self.chevron_height_ratio, 0.01)
        chevron_y = ypos - (chevron_height * self.chevron_vertical_offset_ratio)
        chevron_width = max(self.chevron_min_width_bp, intron_len * self.chevron_width_fraction)
        margin = (chevron_width / 2) + self.chevron_margin_bp
        usable_len = intron_len - (2 * margin)
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

        direction = self._strand_direction(strand)
        start_center = intron_start + margin
        end_center = intron_end - margin
        midpoint = (intron_start + intron_end) / 2.0

        if chevron_count == 1:
            centers = [midpoint]
        else:
            max_symmetric_spacing = (end_center - start_center) / (chevron_count - 1)
            spacing = max_symmetric_spacing
            half_span = spacing * ((chevron_count - 1) / 2.0)
            first_center = midpoint - half_span
            centers = [first_center + (index * spacing) for index in range(chevron_count)]

        for center in centers:
            x_tail = center - (direction * chevron_width / 2)
            x_tip = center + (direction * chevron_width / 2)

            ax.plot(
                [x_tail, x_tip],
                [chevron_y - chevron_height, chevron_y],
                color=self.intron_color,
                linewidth=self.intron_linewidth,
                zorder=1,
            )
            ax.plot(
                [x_tail, x_tip],
                [chevron_y + chevron_height, chevron_y],
                color=self.intron_color,
                linewidth=self.intron_linewidth,
                zorder=1,
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

        strategy = str(self.label_overlap_strategy)
        display_mode = str(self.display)
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
        ax.set_ylim(-0.1, rows * self.row_scale + 0.1)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        self.plot_genes(ax, gr)
        clean_axis(ax)
