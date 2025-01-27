import matplotlib
import pandas as pd
import numpy as np
from coolbox.api import BED, GenomeRange
from coolbox.core.track import Track
from coolbox.core.track.bed.plot import PlotGenes
from coolbox.core.track.bed.base import BedBase
from coolbox.core.track.bed.fetch import FetchBed
from coolbox.utilities.bed import build_bed_index
from loguru import logger as log
import pathlib


class PlotGenesPlotnado(PlotGenes):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.label_location = self.properties.get("label_loc", "mid")
        self.label_y_offset = self.properties.get("label_y_offset", 0)

    def __set_plot_params(self, gr: GenomeRange, ov_genes: pd.DataFrame):
        properties = self.properties
        # bed_type
        self.properties["bed_type"] = properties["bed_type"] or self.infer_bed_type(
            ov_genes
        )
        self.set_colormap(ov_genes)
        # turn labels off when too many intervals are visible.
        if properties["labels"] == "auto":
            if len(ov_genes) > 60:
                self.is_draw_labels = False
            else:
                self.is_draw_labels = True
        self.small_relative = 0.004 * (gr.end - gr.start)
        self.counter = 0

    def __get_length_w(self, fig_width, region_start, region_end):
        """
        to improve the visualization of the genes it is good to have an estimation of the label
        length. In the following code I try to get the length of a 'W' in base pairs.
        """
        if self.is_draw_labels:
            # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
            inches_per_pt = 1.0 / 72.27
            font_in_inches = self.properties["fontsize"] * inches_per_pt
            region_len = region_end - region_start
            bp_per_inch = region_len / fig_width
            font_in_bp = font_in_inches * bp_per_inch
            self.len_w = font_in_bp
            log.debug("len of w set to: {} bp".format(self.len_w))
        else:
            self.len_w = 1

        return self.len_w

    def get_y_pos(self, free_row):
        """
        The y_pos is set such that regions to be plotted do not overlap (stacked). To override this
        the properties['collapsed'] needs to be set.

        The algorithm uses a interval tree (self.region_interval) to check the overlaps
        and a sort of coverage vector 'rows used' to identify the row in which to plot

        Return
        ------
        ypos : int
            y position.
        """

        # if the domain directive is given, ypos simply oscilates between 0 and 100
        if self.properties["display"] == "interlaced":
            return self.properties["interval_height"] if self.counter % 2 == 0 else 1
        elif self.properties["display"] == "collapsed":
            return 0
        else:
            return free_row * self.row_scale

    def plot_genes(
        self, ax, gr: GenomeRange, ov_genes: pd.DataFrame, dry_run=False, fig_width=None
    ):
        properties = self.properties
        self.__set_plot_params(gr, ov_genes)

        assert (not dry_run) or (fig_width is not None)
        if dry_run:
            self.__get_length_w(fig_width, gr.start, gr.end)
        else:
            self.__get_length_w(ax.get_figure().get_figwidth(), gr.start, gr.end)

        num_rows = properties["num_rows"]
        max_num_row_local = 1
        max_ypos = 0
        # check for the number of other intervals that overlap
        #    with the given interval
        #            1         2
        #  012345678901234567890123456
        #  1=========       4=========
        #       2=========
        #         3============
        #
        # for 1 row_last_position = [9]
        # for 2 row_last_position = [9, 14]
        # for 3 row_last_position = [9, 14, 19]
        # for 4 row_last_position = [26, 14, 19]

        row_last_position = []  # each entry in this list contains the end position
        # of genomic interval. The list index is the row
        # in which the genomic interval was plotted.
        # Any new genomic interval that wants to be plotted,
        # knows the row to use by finding the list index that
        # is larger than its start

        # check for overlapping genes including
        # label size (if plotted)

        for bed in ov_genes.itertuples():
            """
            BED12 gene format with exon locations at the end
            chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

            BED9
            bed with rgb at end
            chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44

            BED6
            bed without rgb
            chr2L   0       70000   ID_5    0.26864549832   .

            BED3
            bed with only intervals
            chr2L  0        70000
            """
            self.counter += 1

            if self.is_draw_labels and self.label_location == "right":
                num_name_characters = (
                    len(bed.name) + 2
                )  # +2 to account for an space before and after the name
                bed_extended_end = int(bed.end + (num_name_characters * self.len_w))
            else:
                bed_extended_end = bed.end + 2 * self.small_relative

            # get smallest free row
            if not row_last_position:
                free_row = 0
                row_last_position.append(bed_extended_end)
            else:
                # get list of rows that are less than bed.start, then take the min
                idx_list = [
                    idx
                    for idx, value in enumerate(row_last_position)
                    if value < bed.start
                ]
                if len(idx_list):
                    free_row = min(idx_list)
                    row_last_position[free_row] = bed_extended_end
                else:
                    free_row = len(row_last_position)
                    row_last_position.append(bed_extended_end)

            rgb, edgecolor = self.get_rgb_and_edge_color(bed)
            ypos = self.get_y_pos(free_row)

            # do not plot if the maximum interval rows to plot is reached
            if num_rows and free_row >= float(num_rows):
                continue

            if free_row > max_num_row_local:
                max_num_row_local = free_row
            if ypos > max_ypos:
                max_ypos = ypos

            if not dry_run:
                if properties["bed_type"] == "bed12":
                    if properties["gene_style"] == "flybase":
                        self.draw_gene_with_introns_flybase_style(
                            ax, bed, ypos, rgb, edgecolor
                        )
                    else:
                        self._draw_gene_with_introns(
                            ax,
                            bed,
                            ypos,
                            rgb,
                            edgecolor,
                            arrow_color=properties.get("arrow_color", "black"),
                        )
                else:
                    self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

                self.draw_label(bed, gr, ax, ypos)

        if self.counter == 0:
            log.debug(
                f"*Warning* No intervals were found for file {properties['file']} "
                f"in Track '{properties['name']}' for the interval plotted ({gr}).\n"
            )

        ymax = 0
        if num_rows:
            ymin = float(num_rows) * self.row_scale
            self.current_row_num = num_rows
        else:
            ymin = max_ypos + properties["interval_height"]
            self.current_row_num = len(row_last_position)

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        if not dry_run:
            ax.set_ylim(ymin, ymax)

            if properties["display"] == "collapsed":
                ax.set_ylim(-5, 105)

            ax.set_xlim(gr.start, gr.end)

    def draw_label(self, bed, gr, ax, ypos):
        if self.is_draw_labels:
            # if the label is to be plotted on the right side
            if self.label_location == "right":
                # Check if the label is within the genomic range if not writ the label at the middle of the gene
                if bed.start > gr.start and bed.end < gr.end:
                    self._draw_label_right(bed, gr, ax, ypos)
                elif bed.start < gr.end and bed.end > gr.end:
                    self._draw_label_mid(bed, gr, ax, ypos)

            elif self.label_location == "mid":
                self._draw_label_mid(bed, gr, ax, ypos)

    def _draw_label_right(self, bed, gr, ax, ypos):
        ax.text(
            bed.end + self.small_relative,
            ypos + (float(self.properties["interval_height"]) / 2),
            bed.name,
            horizontalalignment="left",
            verticalalignment="center",
            fontproperties=self.fp,
        )

    def _draw_label_mid(self, bed, gr, ax, ypos):
        # Intersect the gene with the genomic range
        start = max(bed.start, gr.start)
        end = min(bed.end, gr.end)
        # Calculate the middle of the gene
        mid = (start + end) / 2
        ax.text(
            mid,
            ypos
            + (float(self.properties["interval_height"]) / 2)
            + self.label_y_offset,
            bed.name,
            horizontalalignment="center",
            verticalalignment="center",
            fontproperties=self.fp,
        )

    def _draw_gene_with_introns(
        self, ax, bed, ypos, rgb, edgecolor, arrow_color="blue"
    ):
        """
        draws a gene like in flybase gbrowse.
        """
        from matplotlib.patches import Polygon

        properties = self.properties
        height = float(properties["interval_height"])

        if (
            bed.block_count == 0
            and bed.thick_start == bed.start
            and bed.thick_end == bed.end
        ):
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = height / 2
        quarter_height = height / 4
        three_quarter_height = quarter_height * 3

        # draw 'backbone', a line from the start until the end of the gene
        ax.plot(
            [bed.start, bed.end],
            [ypos + half_height, ypos + half_height],
            "black",
            linewidth=0.5,
            zorder=-1,
        )

        for idx in range(bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x1 < bed.thick_start or x0 > bed.thick_end:
                y0 = ypos + quarter_height
                y1 = ypos + three_quarter_height
            else:
                y0 = ypos
                y1 = ypos + height

            if x0 < bed.thick_start < x1:
                vertices = [
                    (x0, ypos + quarter_height),
                    (x0, ypos + three_quarter_height),
                    (bed.thick_start, ypos + three_quarter_height),
                    (bed.thick_start, ypos + height),
                    (bed.thick_start, ypos + height),
                    (x1, ypos + height),
                    (x1, ypos),
                    (bed.thick_start, ypos),
                    (bed.thick_start, ypos + quarter_height),
                ]

            elif x0 < bed.thick_end < x1:
                vertices = [
                    (x0, ypos),
                    (x0, ypos + height),
                    (bed.thick_end, ypos + height),
                    (bed.thick_end, ypos + three_quarter_height),
                    (x1, ypos + three_quarter_height),
                    (x1, ypos + quarter_height),
                    (bed.thick_end, ypos + quarter_height),
                    (bed.thick_end, ypos),
                ]
            else:
                vertices = [(x0, y0), (x0, y1), (x1, y1), (x1, y0)]

            ax.add_patch(
                Polygon(
                    vertices,
                    closed=True,
                    fill=True,
                    linewidth=0.1,
                    edgecolor="none",
                    facecolor=rgb,
                )
            )

            if idx < bed.block_count - 1:
                # plot small arrows using the character '<' or '>' over the back bone
                intron_length = bed.block_starts[idx + 1] - (
                    bed.block_starts[idx] + bed.block_sizes[idx]
                )
                marker = 5 if bed.strand == "+" else 4
                if intron_length > 3 * self.small_relative:
                    pos = np.arange(
                        x1 + 1 * self.small_relative,
                        x1 + intron_length + self.small_relative,
                        int(2 * self.small_relative),
                    )
                    ax.plot(
                        pos,
                        np.zeros(len(pos)) + ypos + half_height,
                        ".",
                        marker=marker,
                        fillstyle="none",
                        color=arrow_color,
                        markersize=3,
                    )

                elif intron_length > self.small_relative:
                    intron_center = x1 + int(intron_length) / 2
                    ax.plot(
                        [intron_center],
                        [ypos + half_height],
                        ".",
                        marker=5,
                        fillstyle="none",
                        color=arrow_color,
                        markersize=3,
                    )


class Genes(BedBase, PlotGenesPlotnado, FetchBed):
    DEFAULT_PROPERTIES = {
        "labels": "on",
    }

    def __init__(
        self,
        file: str = None,
        genome: str = None,
        ignore_file_validation: bool = True,
        min_gene_length: int = 1e4,
        **kwargs,
    ):
        if file is None and genome is None:
            raise ValueError(
                "Genes track requires a file path to a bed file or a genome name"
            )
        elif file is None and genome:
            file = self.get_genes_file(genome)
        elif file:
            file = self.get_genes_file(file)

        properties = BED.DEFAULT_PROPERTIES.copy()
        properties.update(
            {
                "file": file,
                "ignore_file_validation": ignore_file_validation,
                "min_gene_length": min_gene_length,
                **kwargs,
            }
        )
        super().__init__(**properties)
        PlotGenesPlotnado.__init__(self)

        self.min_gene_length = min_gene_length

    def _fetch_intervals(self, file: str, gr: GenomeRange):
        from pybedtools import BedTool

        bt = BedTool(file)
        try:
            bt_tabix = bt.tabix(force=True)
            intervals = bt_tabix.tabix_intervals(f"{gr.chrom}:{gr.start}-{gr.end}")

        except OSError:  # Handle the case where the bed file is not tabix indexed or the user does not have permission to write to the directory
            import tempfile

            with tempfile.NamedTemporaryFile() as tmp:
                bt.saveas(tmp.name)
                bt_tabix = BedTool(tmp.name).tabix(force=True)
                intervals = bt_tabix.tabix_intervals(f"{gr.chrom}:{gr.start}-{gr.end}")

        if len(intervals) == 0:
            log.warning(
                f"No intervals were found for file {file} in the interval plotted ({gr}).\n"
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

        df['block_starts'] = df['block_starts'].apply(lambda x: list(map(int, x.split(','))))
        df['block_sizes'] = df['block_sizes'].apply(lambda x: list(map(int, x.split(','))))

        return df
                                                    

    def fetch_data(self, gr: GenomeRange, **kwargs):
        intervals = self._fetch_intervals(self.properties["file"], gr)
        intervals = intervals.query("(end - start) > @self.min_gene_length")
        return intervals

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        ov_intervals: pd.DataFrame = self.fetch_plot_data(gr, **kwargs)
        self.plot_genes(ax, gr, ov_intervals)

    def get_genes_file(self, file: str):
        import importlib
        import json

        try:
            bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
            bed_paths = bed_prefix / "genes.json"

            with open(bed_paths) as f:
                gene_files = json.load(f)

        except FileNotFoundError:
            if pathlib.Path(file).exists():
                return file
            else:
                raise FileNotFoundError(f"File {file} does not exist")

        if file in gene_files:
            return bed_prefix / gene_files[file]
        else:
            return file
