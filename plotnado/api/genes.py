import matplotlib
import pandas as pd
from coolbox.api import BED, GenomeRange
from coolbox.core.track import Track
from coolbox.core.track.bed.base import BedBase
from coolbox.core.track.bed.fetch import FetchBed
from coolbox.core.track.bed.plot import PlotBed
from coolbox.utilities.bed import build_bed_index
from loguru import logger as log


class PlotGenes(PlotBed):
    def plot_genes(self, ax, gr: GenomeRange, ov_genes: pd.DataFrame):
        properties = self.properties
        # bed_type
        self.properties["bed_type"] = properties["bed_type"] or self.infer_bed_type(
            ov_genes
        )
        # as min_score and max_score change every plot, we compute them for every plot
        min_score, max_score = properties["min_score"], properties["max_score"]
        has_score_col = properties["bed_type"] in ("bed6", "bed9", "bed12")
        if has_score_col and len(ov_genes):
            min_score = (min_score != "inf") or ov_genes["score"].min()
            max_score = (max_score != "-inf") or ov_genes["score"].max()
        min_score, max_score = float(min_score), float(max_score)

        # set colormap
        if self.colormap is not None:
            norm = matplotlib.colors.Normalize(vmin=min_score, vmax=max_score)
            cmap = matplotlib.cm.get_cmap(properties["color"])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        if properties["color"] == "bed_rgb" and properties["bed_type"] not in [
            "bed12",
            "bed9",
        ]:
            log.warning(
                "*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                "been set to {}".format(self.COLOR)
            )
            self.properties["color"] = self.COLOR
            self.colormap = None

        self.counter = 0
        self.small_relative = 0.004 * (gr.end - gr.start)
        self.get_length_w(ax.get_figure().get_figwidth(), gr.start, gr.end)
        # turn labels off when too many intervals are visible.
        if properties["labels"] == "on" and len(ov_genes) > 60:
            self.is_draw_labels = False

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

            if self.is_draw_labels:
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

            if properties["bed_type"] == "bed12":
                if properties["gene_style"] == "flybase":
                    self.draw_gene_with_introns_flybase_style(
                        ax, bed, ypos, rgb, edgecolor
                    )
                else:
                    self.draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor)
            else:
                self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

            if (
                self.is_draw_labels
                and (bed.start >= gr.start and bed.start <= gr.end)
                or (bed.end >= gr.start and bed.end <= gr.end)
            ):
                overlap_start = max(bed.start, gr.start)
                overlap_end = min(bed.end, gr.end)
                overlap_center = (overlap_start + overlap_end) / 2

                ax.text(
                    overlap_center + self.small_relative,
                    ypos - (float(properties["interval_height"]) / 4),
                    bed.name,
                    horizontalalignment="left",
                    verticalalignment="center",
                    fontproperties=self.fp,
                )

        if self.counter == 0:
            log.warning(
                f"*Warning* No intervals were found for file {properties['file']} "
                f"in Track '{properties['name']}' for the interval plotted ({gr}).\n"
            )

        ymax = 0
        if num_rows:
            ymin = float(num_rows) * self.row_scale
        else:
            ymin = max_ypos + properties["interval_height"]

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

        if properties["display"] == "domain":
            ax.set_ylim(-5, 205)
        elif properties["display"] == "collapsed":
            ax.set_ylim(-5, 105)

        ax.set_xlim(gr.start, gr.end)


class GenesBase(Track, PlotGenes):
    """
    BED Base track.

    Parameters
    ----------
    style : {'gene', 'tad'}

    gene_style: {'flybase', 'normal'}

    display : {'stacked', 'interlaced', 'collapsed'}, optional
        Display mode. (Default: 'stacked')

    color : str, optional
        Track color, 'bed_rgb' for auto specify color according to bed record.
        (Default: 'bed_rgb')

    border_color : str, optional
        Border_color of gene. (Default: 'black')

    fontsize : int, optional
        Font size. (Default: BED.DEFAULT_FONTSIZE)

    labels : {True, False, 'auto'}, optional
        Draw bed name or not. 'auto' for automate decision according to density.
        (Default: 'auto')

    interval_height : int, optional
        The height of the interval. (Default: 100)

    num_rows : int, optional
        Set the max interval rows. (Default: unlimited interval rows)

    max_value : float, optional
        Max score. (Default: inf)

    min_value : float, optional
        Min score. (Default: -inf)

    border_style: str, optional
        Border style of tad. (Default: 'solid')

    border_width: int, optional
        Border width of tad. (Default: '2.0')

    show_score : bool
        Show bed score or not.
        default False.

    score_font_size : {'auto', int}
        Score text font size.
        default 'auto'

    score_font_color : str
        Score text color.
        default '#000000'

    score_height_ratio : float
        (text tag height) / (TAD height). used for adjust the position of Score text.
        default 0.5

    border_only : bool
        Only show border, default False

    """

    STYLE_GENE = "gene"
    STYLE_TAD = "tad"

    COLOR = "#1f78b4"

    DEFAULT_PROPERTIES = {
        "style": STYLE_GENE,
        # gene
        "gene_style": "flybase",
        "display": "stacked",
        "color": "bed_rgb",
        "border_color": "#1f78b4",
        "fontsize": 12,
        "interval_height": 100,
        "num_rows": None,
        "labels": "off",
        "min_score": "-inf",
        "max_score": "inf",
        "bed_type": None,
        # tad
        "border_style": "--",
        "border_width": 2.0,
        "show_score": False,
        "score_font_size": "auto",
        "score_font_color": "#000000",
        "score_height_ratio": 0.4,
        "border_only": False,
    }

    def __init__(self, **kwargs):
        properties = BedBase.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(properties)
        self.init_for_plot()

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            BED interval table. The table should be in format like:

            bed_fields = ['chromosome', 'start', 'end',
                          'name', 'score', 'strand',
                          'thick_start', 'thick_end',
                          'rgb', 'block_count',
                          'block_sizes', 'block_starts']
            The table can be in bed6/bed9/bed12 format and the trailing columns can be omited.

        """
        raise NotImplementedError

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        ov_intervals: pd.DataFrame = self.fetch_plot_data(gr, **kwargs)

        style = self.properties["style"]
        if style == self.STYLE_TAD:
            self.plot_tads(ax, gr, ov_intervals)
        elif style == self.STYLE_GENE:
            self.plot_genes(ax, gr, ov_intervals)
        else:
            raise ValueError("style not supportted, should be one of 'gene' 'tad' ")
        self.plot_label()


class Genes(GenesBase, FetchBed):
    """
    Bed Track for plotting 1d intervals data from .bed file.
    The input bed file can be bed3/bed6/bed9/bed12

    Parameters
    ----------
    file: str
        The file path of `.bed` file.


    """

    DEFAULT_PROPERTIES = {
        "labels": "on",
    }

    def __init__(self, file, **kwargs):
        
        file = self.get_genes_file(file)

        properties = BED.DEFAULT_PROPERTIES.copy()
        properties.update({"file": file, **kwargs})
        super().__init__(**properties)
        self.bgz_file = build_bed_index(file)

    def fetch_data(self, gr: GenomeRange, **kwargs):
        return self.fetch_intervals(self.bgz_file, gr)

    def get_genes_file(self, file: str):

        import importlib
        import json

        bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
        bed_paths = bed_prefix / "genes.json"

        with open(bed_paths) as f:
            gene_files = json.load(f)

        if file in gene_files:
            return bed_prefix / gene_files[file]
        else:
            return file 
