"""
Genes track for displaying gene annotations.
"""

import importlib.resources
import json
from pathlib import Path
from typing import Literal

import matplotlib.axes
import matplotlib.patches
import pandas as pd
from pydantic import BaseModel, ConfigDict

from .base import Track
from .region import GenomicRegion
from .utils import clean_axis, read_bed_regions, read_gtf_regions


class GenesAesthetics(BaseModel):
    """Aesthetics configuration for gene tracks."""

    style: Literal["std"] = "std"
    color: str = "#2c3e50"
    fill: bool = True
    alpha: float = 1.0
    display: Literal["collapsed", "expanded"] = "collapsed"
    minimum_gene_length: int = 0
    max_number_of_rows: int = 4
    interval_height: float = 0.1
    arrow_size: float = 0.12
    exon_linewidth: float = 0.8
    exon_edge_color: str = "black"
    intron_linewidth: float = 0.8
    intron_color: str = "black"
    gene_label_font_size: int = 6
    gene_label_style: Literal["normal", "italic", "oblique"] = "italic"


class Genes(Track):
    """Track for displaying gene annotations from BED12 or GTF files."""

    data: Path | str | pd.DataFrame | None = None
    genome: str | None = None
    aesthetics: GenesAesthetics = GenesAesthetics()
    height: float = 1.5

    row_scale: float = 1.0
    small_relative: float = 0.01
    gene_count: int = 0

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

        rename_map = {
            "name": "geneid",
            "Name": "geneid",
            "itemRgb": "rgb",
            "blockCount": "block_count",
            "blockSizes": "block_sizes",
            "blockStarts": "block_starts",
            "thickStart": "thick_start",
            "thickEnd": "thick_end",
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        if "geneid" not in df.columns:
            df["geneid"] = "gene"
        if "strand" not in df.columns:
            df["strand"] = "+"
        if "block_count" not in df.columns:
            df["block_count"] = 1
            df["block_sizes"] = (df["end"] - df["start"]).astype(int).astype(str)
            df["block_starts"] = "0"

        df["block_starts"] = df["block_starts"].apply(self._parse_int_list)
        df["block_sizes"] = df["block_sizes"].apply(self._parse_int_list)
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
                facecolor=self.color,
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

    def _draw_gene_simple(
        self, ax: matplotlib.axes.Axes, gene: pd.Series, ypos: float
    ) -> None:
        rect = matplotlib.patches.Rectangle(
            (gene.start, ypos - self.interval_height / 2),
            gene.end - gene.start,
            self.interval_height,
            linewidth=self.exon_linewidth,
            edgecolor=self.exon_edge_color,
            facecolor=self.color,
            alpha=self.alpha,
        )
        ax.add_patch(rect)

    def _draw_gene_feature(
        self, ax: matplotlib.axes.Axes, gene: pd.Series, row_index: int
    ) -> None:
        ypos = self.get_y_pos(row_index)
        block_count = int(gene.get("block_count", 1))
        if block_count > 1:
            self._draw_gene_with_introns(ax, gene, ypos)
        else:
            self._draw_gene_simple(ax, gene, ypos)

        label_offset = max(self.small_relative * 2, 1000)
        ax.text(
            gene["end"] + label_offset,
            ypos,
            str(gene.get("geneid", "gene")),
            va="center",
            ha="left",
            fontsize=self.gene_label_font_size,
            style=self.gene_label_style,
            color=self.color,
        )

    def plot_genes(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        genes_df = self.fetch_data(gr)

        if genes_df.empty:
            ax.set_ylim(0, 1)
            ax.set_xlim(gr.start, gr.end)
            return

        row_last_positions: list[int] = []
        max_row_index = 0

        if self.display == "expanded":
            for _, gene in genes_df.iterrows():
                row_index = self._allocate_row_index(
                    row_last_positions, int(gene["start"]), int(gene["end"])
                )
                row_index = min(row_index, self.max_number_of_rows - 1)
                max_row_index = max(max_row_index, row_index)
                self._draw_gene_feature(ax, gene, row_index)
        else:
            for _, gene in genes_df.iterrows():
                self._draw_gene_feature(ax, gene, 0)

        ax.set_xlim(gr.start, gr.end)
        rows = max(1, max_row_index + 1)
        ax.set_ylim(-0.1, rows * self.row_scale + 0.1)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        self.plot_genes(ax, gr)
        clean_axis(ax)
