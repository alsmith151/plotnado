"""
Genes track for displaying gene annotations.
"""

import importlib.resources
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Literal, Optional, Union

import matplotlib.axes
import matplotlib.patches
import pandas as pd
from pydantic import BaseModel

from .region import GenomicRegion
from .base import Track
from .utils import clean_axis


class GenesAesthetics(BaseModel):
    """
    Aesthetics configuration for gene tracks.
    """

    style: Literal["std"] = "std"
    color: str = "black"
    fill: bool = True
    alpha: float = 1.0
    display: Literal["collapsed", "expanded"] = "collapsed"
    minimum_gene_length: int = 0
    max_number_of_rows: int = 4
    interval_height: float = 0.5
    arrow_size: float = 0.1
    title: str = "Genes"


def tabix_gtf(file: Path) -> Path:
    """Create a tabix index for a GTF file."""
    sorted_file = file.parent / f"{file.name}.sorted"
    gz_file = file.parent / f"{file.name}.gz"

    with sorted_file.open("w") as sf:
        subprocess.run(
            ["sort", "-k1,1", "-k4,4n", str(file)],
            stdout=sf,
            check=True,
        )

    with gz_file.open("wb") as gf:
        subprocess.run(
            ["bgzip", "-c", str(sorted_file)],
            stdout=gf,
            check=True,
        )

    subprocess.run(["tabix", "-p", "gff", str(gz_file)], check=True)
    sorted_file.unlink()
    return gz_file


def gtf_line_to_bed12_line(df: pd.DataFrame) -> pd.Series:
    """Convert GTF records to BED12 format."""
    df = df.sort_values(["seqname", "start"])
    geneid = df["geneid"].iloc[0]
    exons = df.query('feature == "exon"')
    chrom = df["seqname"].iloc[0]
    start = str(df["start"].min())
    end = str(df["end"].max())
    strand = df["strand"].iloc[0]
    thick_start = start if strand == "+" else end
    thick_end = thick_start
    color = "0,0,0"
    block_count = str(exons.shape[0])
    block_sizes = ",".join((exons["end"] - exons["start"]).values.astype(str))
    block_starts = ",".join((exons["start"] - int(start)).astype(str))

    return pd.Series(
        {
            "chrom": chrom,
            "start": start,
            "end": end,
            "geneid": geneid,
            "score": "0",
            "strand": strand,
            "thick_start": thick_start,
            "thick_end": thick_end,
            "color": color,
            "block_count": block_count,
            "block_sizes": block_sizes,
            "block_starts": block_starts,
        }
    )


class Genes(Track):
    """
    Track for displaying gene annotations from BED12 or GTF files.

    Attributes:
        data: Path to BED12 or GTF file (optional if genome is provided)
        genome: Genome name (e.g., 'hg38') to use built-in annotations
        aesthetics: Visual styling configuration
    """

    data: Optional[Union[Path, str, pd.DataFrame]] = None
    genome: Optional[str] = None
    aesthetics: GenesAesthetics = GenesAesthetics()
    height: float = 1.5

    row_scale: float = 0.2
    small_relative: float = 0.01
    len_w: float = 0.01
    gene_count: int = 0
    current_row_num: int = 0

    class Config:
        arbitrary_types_allowed = True

    def _fetch_genes_from_package(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch genes from built-in package data."""
        try:
            bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
            bed_paths = bed_prefix / "genes.json"

            with open(bed_paths) as f:
                gene_files = json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(
                "The genes database is not available. Please provide a data file."
            )

        if self.genome not in gene_files:
            raise ValueError(
                f"Genome {self.genome} not found in the genes database. "
                f"Available genomes: {list(gene_files.keys())}"
            )

        gene_file = bed_prefix / gene_files[self.genome]
        return self._fetch_from_disk_bed12(gr, gene_file)

    def _fetch_from_disk_bed12(
        self, gr: GenomicRegion, file_path: Optional[Path] = None
    ) -> pd.DataFrame:
        """Fetch genes from a BED12 file."""
        from pybedtools import BedTool

        path = str(file_path or self.data)
        bt = BedTool(path)

        try:
            bt_tabix = bt.tabix(force=True)
            intervals = bt_tabix.tabix_intervals(f"{gr.chromosome}:{gr.start}-{gr.end}")
        except OSError:
            with tempfile.NamedTemporaryFile() as tmp:
                bt.saveas(tmp.name)
                bt_tabix = BedTool(tmp.name).tabix(force=True)
                intervals = bt_tabix.tabix_intervals(
                    f"{gr.chromosome}:{gr.start}-{gr.end}"
                )

        df = intervals.to_dataframe()
        df = df.rename(
            columns={
                "itemRgb": "rgb",
                "blockCount": "block_count",
                "blockSizes": "block_sizes",
                "blockStarts": "block_starts",
                "thickStart": "thick_start",
                "thickEnd": "thick_end",
            }
        )

        if not df.empty:
            df["block_starts"] = df["block_starts"].apply(
                lambda x: list(map(int, x.split(",")))
            )
            df["block_sizes"] = df["block_sizes"].apply(
                lambda x: list(map(int, x.split(",")))
            )

        return df

    def _fetch_from_disk_gtf(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch genes from a GTF file."""
        import pysam

        data_path = str(self.data)

        if (
            not data_path.endswith(".gz")
            and not Path(data_path).with_suffix(".gz.tbi").exists()
        ):
            try:
                import plotnado.tracks as pg_tracks

                data = pg_tracks.tabix_gtf(Path(data_path))
            except OSError:
                temp_file = tempfile.NamedTemporaryFile(delete=False)
                shutil.copy(data_path, temp_file.name)
                import plotnado.tracks as pg_tracks

                data = pg_tracks.tabix_gtf(Path(temp_file.name))
        else:
            data = data_path

        gtf = pysam.TabixFile(data)
        records = []
        for record in gtf.fetch(gr.chromosome, gr.start, gr.end, parser=pysam.asGTF()):
            # Handle both pysam objects and dictionaries (for mocks)
            rec_dict = {}
            for field in [
                "seqname",
                "start",
                "end",
                "strand",
                "feature",
                "attributes",
                "attribute",
            ]:
                if hasattr(record, field):
                    rec_dict[field] = getattr(record, field)
                elif isinstance(record, dict) and field in record:
                    rec_dict[field] = record[field]
            records.append(rec_dict)

        df = pd.DataFrame(records)
        df = df.rename(columns={"seqname": "chrom"})
        if df.empty:
            return pd.DataFrame(
                columns=["chrom", "start", "end", "strand", "feature", "geneid"]
            )

        # Consolidate attribute/attributes
        if "attributes" in df.columns:
            df["attribute"] = df["attributes"]
            df = df.drop(columns=["attributes"])

        # If 'attribute' is still duplicated (as a column name), take the first one
        if isinstance(df.get("attribute"), pd.DataFrame):
            df["attribute"] = df["attribute"].iloc[:, 0]

        if "attribute" in df.columns:
            df["geneid"] = (
                df["attribute"].astype(str).str.extract(r'gene_id\s?"(.*?)";.*')
            )
        else:
            df["geneid"] = None
        df = df.query('feature.isin(["5UTR", "3UTR", "exon"])')
        df = df.loc[lambda df: df["chrom"].str.contains(r"^chr[xXYy]?[1-9]?[0-9]?$")]

        records = []
        for gene, gene_df in df.sort_values(["chrom", "start"]).groupby("geneid"):
            bed12_line = gtf_line_to_bed12_line(gene_df)
            records.append(bed12_line)

        return pd.DataFrame(records)

    def fetch_data(self, gr: GenomicRegion) -> pd.DataFrame:
        """Fetch gene data for the given region."""
        if self.data is None and self.genome is None:
            raise ValueError("Either data or genome must be provided")
        elif self.data is None:
            data = self._fetch_genes_from_package(gr)
        elif isinstance(self.data, pd.DataFrame):
            data = self.data
        else:
            fn_string = str(self.data)
            if fn_string.endswith(".gtf.gz") or fn_string.endswith(".gtf"):
                data = self._fetch_from_disk_gtf(gr)
            elif fn_string.endswith(".bed.gz") or fn_string.endswith(".bed"):
                data = self._fetch_from_disk_bed12(gr)
            else:
                raise ValueError(
                    "Unsupported file format. Only .gtf and .bed files are supported."
                )

        # Filter by minimum gene length
        if self.aesthetics.minimum_gene_length > 0:
            data = data.query(f"end - start >= {self.aesthetics.minimum_gene_length}")
        return data

    def _compute_extended_end_bp(
        self, end_bp: int, label: str, char_width_bp: float
    ) -> int:
        """Compute gene end position extended by label padding."""
        label_len = len(label)
        padding_bp = label_len * char_width_bp + 2 * self.small_relative
        return int(end_bp + padding_bp)

    def _estimate_label_char_width(self, ax: matplotlib.axes.Axes) -> float:
        """Estimate the width in bp of a single character for gene labels."""
        # Simple estimation based on figure width and axis range
        fig_width_inches = ax.get_figure().get_size_inches()[0]
        ax_range_bp = ax.get_xlim()[1] - ax.get_xlim()[0]
        # Rough estimate: 1 char is 0.1 inches at default font size
        char_width_inches = 0.1
        return (char_width_inches / fig_width_inches) * ax_range_bp

    def _draw_gene_feature(
        self, ax: matplotlib.axes.Axes, gene: pd.Series, row_index: int
    ) -> None:
        """Draw a single gene feature on the axes."""
        ypos = self.get_y_pos(row_index)
        fill_color = self.aesthetics.color
        edge_color = "black"

        # Check if we should draw with introns
        if "block_count" in gene and int(gene["block_count"]) > 1:
            self._draw_gene_with_introns(ax, gene, ypos, fill_color, edge_color)
        else:
            self._draw_gene_simple(ax, gene, ypos, fill_color, edge_color)

        # Add label
        if "geneid" in gene:
            ax.text(
                gene["end"] + self.small_relative,
                ypos,
                str(gene["geneid"]),
                va="center",
                ha="left",
                fontsize=8,
            )

    def _allocate_row_index(
        self, row_last_positions: List[int], start_bp: int, end_bp: int
    ) -> int:
        """Allocate a row index where the new feature does not overlap existing ones."""
        for idx, last_end in enumerate(row_last_positions):
            if last_end < start_bp:
                row_last_positions[idx] = end_bp
                return idx
        row_last_positions.append(end_bp)
        return len(row_last_positions) - 1

    def _draw_gene_with_introns(
        self, ax, gene, ypos: float, fill_color, edge_color
    ) -> None:
        """Draw a gene with exon/intron structure."""
        if not hasattr(gene, "block_count") or gene.block_count == 0:
            self._draw_gene_simple(ax, gene, ypos, fill_color, edge_color)
            return

        block_starts = gene.block_starts if isinstance(gene.block_starts, list) else []
        block_sizes = gene.block_sizes if isinstance(gene.block_sizes, list) else []

        # Draw intron lines
        for i in range(len(block_starts) - 1):
            intron_start = gene.start + block_starts[i] + block_sizes[i]
            intron_end = gene.start + block_starts[i + 1]
            ax.plot(
                [intron_start, intron_end], [ypos, ypos], color=edge_color, linewidth=1
            )

        # Draw exons
        for start, size in zip(block_starts, block_sizes):
            exon_start = gene.start + start
            rect = matplotlib.patches.Rectangle(
                (exon_start, ypos - self.aesthetics.interval_height / 2),
                size,
                self.aesthetics.interval_height,
                linewidth=1,
                edgecolor=edge_color,
                facecolor=fill_color,
                alpha=self.aesthetics.alpha,
            )
            ax.add_patch(rect)

        # Draw strand arrow
        if hasattr(gene, "strand"):
            arrow_y = ypos
            if gene.strand == "+":
                arrow_x = gene.end
                ax.arrow(
                    arrow_x - self.aesthetics.arrow_size,
                    arrow_y,
                    self.aesthetics.arrow_size,
                    0,
                    head_width=self.aesthetics.interval_height / 2,
                    head_length=self.aesthetics.arrow_size / 2,
                    fc=edge_color,
                    ec=edge_color,
                )
            else:
                arrow_x = gene.start
                ax.arrow(
                    arrow_x + self.aesthetics.arrow_size,
                    arrow_y,
                    -self.aesthetics.arrow_size,
                    0,
                    head_width=self.aesthetics.interval_height / 2,
                    head_length=self.aesthetics.arrow_size / 2,
                    fc=edge_color,
                    ec=edge_color,
                )

    def _draw_gene_simple(self, ax, gene, ypos: float, fill_color, edge_color) -> None:
        """Draw a simple gene representation."""
        rect = matplotlib.patches.Rectangle(
            (gene.start, ypos - self.aesthetics.interval_height / 2),
            gene.end - gene.start,
            self.aesthetics.interval_height,
            linewidth=1,
            edgecolor=edge_color,
            facecolor=fill_color,
            alpha=self.aesthetics.alpha,
        )
        ax.add_patch(rect)

    def get_y_pos(self, row_index: int) -> float:
        """Get the y position for a given row index."""
        return row_index * self.row_scale

    def plot_genes(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Render gene features along a genomic range."""
        overlapping_genes = self.fetch_data(gr)

        import logging

        logger = logging.getLogger("plotnado")

        if overlapping_genes.empty:
            logger.warning(f"No genes found in {gr}")
            ax.set_ylim(0, 1)
            ax.set_xlim(gr.start, gr.end)
            return

        row_last_positions: List[int] = []
        self.gene_count = 0

        # Estimate row count if expanding
        if self.aesthetics.display == "expanded":
            label_char_width = self._estimate_label_char_width(ax)
            for _, gene in overlapping_genes.iterrows():
                extended_end = self._compute_extended_end_bp(
                    gene["end"], gene["geneid"], label_char_width
                )
                row_index = self._allocate_row_index(
                    gene["start"], extended_end, row_last_positions
                )
                self.gene_count += 1
                self._draw_gene_feature(ax, gene, row_index)
        else:
            # Collapsed mode
            for _, gene in overlapping_genes.iterrows():
                self.gene_count += 1
                self._draw_gene_feature(ax, gene, 0)

        ax.set_xlim(gr.start, gr.end)

    def plot(self, ax: matplotlib.axes.Axes, gr: GenomicRegion) -> None:
        """Plot the Genes track."""
        self.plot_genes(ax, gr)
        clean_axis(ax)
