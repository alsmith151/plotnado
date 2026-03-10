"""Tests for the refactored Genes track."""

from unittest.mock import MagicMock, patch

import matplotlib.markers
import matplotlib.pyplot as plt
import pandas as pd

from plotnado.tracks import Genes, GenesAesthetics, GenomicRegion


def _chevron_calls(ax: MagicMock) -> list:
    return [
        call
        for call in ax.plot.call_args_list
        if call.kwargs.get("marker")
        in {matplotlib.markers.CARETRIGHTBASE, matplotlib.markers.CARETLEFTBASE}
    ]


class TestGenes:
    def test_init(self):
        aesthetics = GenesAesthetics(color="blue", display="collapsed")
        genes = Genes(title="Test Genes", aesthetics=aesthetics)

        assert genes.title == "Test Genes"
        assert genes.aesthetics.color == "blue"
        assert genes.aesthetics.display == "collapsed"

    @patch("plotnado.tracks.genes.read_gtf_regions")
    def test_fetch_from_disk_gtf(self, mock_read_gtf, genomic_region):
        mock_read_gtf.return_value = pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [1100, 1300],
                "End": [1200, 1400],
                "Feature": ["exon", "exon"],
                "gene_id": ["gene1", "gene1"],
                "Strand": ["+", "+"],
            }
        )

        genes = Genes(data="/path/to/mock.gtf")
        result = genes._fetch_from_disk_gtf(genomic_region)

        mock_read_gtf.assert_called_once_with(
            "/path/to/mock.gtf",
            genomic_region.chromosome,
            genomic_region.start,
            genomic_region.end,
        )
        assert isinstance(result, pd.DataFrame)
        assert set(["chrom", "start", "end", "geneid"]).issubset(result.columns)

    @patch("plotnado.tracks.genes.read_bed_regions")
    def test_fetch_from_disk_bed12(self, mock_read_bed, genomic_region):
        mock_read_bed.return_value = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1300],
                "end": [1200, 1400],
                "name": ["gene1", "gene2"],
                "strand": ["+", "+"],
                "blockCount": [1, 1],
                "blockSizes": ["100", "100"],
                "blockStarts": ["0", "0"],
            }
        )

        genes = Genes(data="/path/to/mock.bed")
        result = genes._fetch_from_disk_bed12(genomic_region)

        mock_read_bed.assert_called_once_with(
            "/path/to/mock.bed",
            genomic_region.chromosome,
            genomic_region.start,
            genomic_region.end,
        )
        assert "block_starts" in result.columns
        assert "block_sizes" in result.columns
        assert "geneid" in result.columns

    @patch("plotnado.tracks.genes.read_bed_regions")
    def test_fetch_from_disk_bed12_parses_capitalized_bed12_columns(
        self, mock_read_bed, genomic_region
    ):
        mock_read_bed.return_value = pd.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [1000],
                "End": [1600],
                "Name": ["gene1"],
                "Strand": ["-"],
                "BlockCount": [2],
                "BlockSizes": ["100,150"],
                "BlockStarts": ["0,450"],
            }
        ).rename(columns={"Chromosome": "chrom", "Start": "start", "End": "end"})

        genes = Genes(data="/path/to/mock.bed")
        result = genes._fetch_from_disk_bed12(genomic_region)

        assert int(result.iloc[0]["block_count"]) == 2
        assert result.iloc[0]["strand"] == "-"
        assert result.iloc[0]["block_starts"] == [0, 450]
        assert result.iloc[0]["block_sizes"] == [100, 150]

    @patch("plotnado.tracks.genes.read_bed_regions")
    def test_fetch_from_disk_bed12_parses_bigbed_field_columns(
        self, mock_read_bed, genomic_region
    ):
        mock_read_bed.return_value = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1000],
                "end": [1800],
                "field_1": ["gene1"],
                "field_2": [0],
                "field_3": ["+"],
                "field_4": [1000],
                "field_5": [1800],
                "field_6": ["0"],
                "field_7": [2],
                "field_8": ["100,200"],
                "field_9": ["0,600"],
            }
        )

        genes = Genes(data="/path/to/mock.bigbed")
        result = genes._fetch_from_disk_bed12(genomic_region)

        assert int(result.iloc[0]["block_count"]) == 2
        assert result.iloc[0]["block_starts"] == [0, 600]
        assert result.iloc[0]["block_sizes"] == [100, 200]

    def test_fetch_data_dataframe_filters_length(self, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr1"],
                "start": [1100, 1300, 1500],
                "end": [1150, 1400, 1600],
                "geneid": ["gene1", "gene2", "gene3"],
                "block_count": [1, 1, 1],
                "block_sizes": [[50], [100], [100]],
                "block_starts": [[0], [0], [0]],
            }
        )

        genes = Genes(data=data, aesthetics=GenesAesthetics(minimum_gene_length=80))
        result = genes.fetch_data(genomic_region)

        assert result.shape[0] == 2
        assert set(result["geneid"]) == {"gene2", "gene3"}

    def test_plot_genes_no_data(self, mock_ax, genomic_region):
        genes = Genes(data=pd.DataFrame(columns=["chrom", "start", "end"]))
        genes.plot_genes(mock_ax, genomic_region)

        mock_ax.set_ylim.assert_called_once_with(0, 1)
        mock_ax.set_xlim.assert_called_once_with(genomic_region.start, genomic_region.end)

    @patch.object(Genes, "_draw_gene_feature")
    def test_plot_genes_expanded(self, mock_draw, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1150],
                "end": [1300, 1450],
                "geneid": ["gene1", "gene2"],
                "strand": ["+", "+"],
                "block_count": [1, 1],
                "block_sizes": [[200], [300]],
                "block_starts": [[0], [0]],
            }
        )

        genes = Genes(data=data, aesthetics=GenesAesthetics(display="expanded"))
        genes.plot_genes(mock_ax, genomic_region)

        assert mock_draw.call_count == 2

    def test_fetch_data_requires_source(self, genomic_region):
        genes = Genes()
        try:
            genes.fetch_data(genomic_region)
            assert False, "Expected ValueError"
        except ValueError as exc:
            assert "Either data or genome must be provided" in str(exc)

    def test_draw_gene_with_introns_adds_chevrons(self):
        ax = MagicMock()
        genes = Genes(data=pd.DataFrame())
        gene = pd.Series(
            {
                "start": 1000,
                "end": 6000,
                "strand": "+",
                "block_starts": [0, 4000],
                "block_sizes": [400, 500],
                "block_count": 2,
            }
        )

        genes._draw_gene_with_introns(ax, gene, ypos=0.5)

        assert len(_chevron_calls(ax)) >= 1

    def test_draw_gene_with_introns_chevrons_follow_minus_strand(self):
        ax = MagicMock()
        genes = Genes(data=pd.DataFrame())
        gene = pd.Series(
            {
                "start": 1000,
                "end": 6000,
                "strand": "-1",
                "block_starts": [0, 4000],
                "block_sizes": [400, 500],
                "block_count": 2,
            }
        )

        genes._draw_gene_with_introns(ax, gene, ypos=0.5)

        assert _chevron_calls(ax)[0].kwargs["marker"] == matplotlib.markers.CARETLEFTBASE

    def test_draw_gene_with_introns_chevrons_centered_on_intron_midpoint(self):
        ax = MagicMock()
        genes = Genes(data=pd.DataFrame())
        gene = pd.Series(
            {
                "start": 1000,
                "end": 2600,
                "strand": "+",
                "block_starts": [0, 1200],
                "block_sizes": [300, 300],
                "block_count": 2,
            }
        )

        genes._draw_gene_with_introns(ax, gene, ypos=0.5)

        intron_start = 1000 + 300
        intron_end = 1000 + 1200
        intron_mid = (intron_start + intron_end) / 2
        first_chevron_x = _chevron_calls(ax)[0].args[0]
        chevron_mid = first_chevron_x[0]

        assert abs(chevron_mid - intron_mid) < 1e-6

    def test_draw_gene_with_introns_chevrons_evenly_spaced_on_long_intron(self):
        ax = MagicMock()
        genes = Genes(data=pd.DataFrame())
        gene = pd.Series(
            {
                "start": 1000,
                "end": 52000,
                "strand": "+",
                "block_starts": [0, 50000],
                "block_sizes": [500, 500],
                "block_count": 2,
            }
        )

        genes._draw_gene_with_introns(ax, gene, ypos=0.5)

        chevron_line_calls = _chevron_calls(ax)
        centers = [call.args[0][0] for call in chevron_line_calls]

        assert len(centers) >= 6
        spacings = [centers[index + 1] - centers[index] for index in range(len(centers) - 1)]
        assert spacings
        assert max(spacings) - min(spacings) < 1e-6

    def test_draw_gene_with_introns_has_center_chevron_at_intron_midpoint(self):
        ax = MagicMock()
        genes = Genes(data=pd.DataFrame())
        gene = pd.Series(
            {
                "start": 1000,
                "end": 91000,
                "strand": "+",
                "block_starts": [0, 90000],
                "block_sizes": [500, 500],
                "block_count": 2,
            }
        )

        genes._draw_gene_with_introns(ax, gene, ypos=0.5)

        intron_start = 1000 + 500
        intron_end = 1000 + 90000
        intron_mid = (intron_start + intron_end) / 2
        chevron_line_calls = _chevron_calls(ax)
        centers = [call.args[0][0] for call in chevron_line_calls]

        assert any(abs(center - intron_mid) < 1e-6 for center in centers)

    def test_draw_gene_with_introns_respects_chevron_target_spacing(self):
        ax = MagicMock()
        genes = Genes(
            data=pd.DataFrame(),
            aesthetics=GenesAesthetics(chevron_target_spacing_bp=20_000),
        )
        gene = pd.Series(
            {
                "start": 1000,
                "end": 91000,
                "strand": "+",
                "block_starts": [0, 90000],
                "block_sizes": [500, 500],
                "block_count": 2,
            }
        )

        genes._draw_gene_with_introns(ax, gene, ypos=0.5)

        chevron_line_calls = _chevron_calls(ax)
        assert len(chevron_line_calls) < 6

    def test_draw_gene_with_introns_respects_chevron_vertical_offset(self):
        ax = MagicMock()
        genes = Genes(
            data=pd.DataFrame(),
            aesthetics=GenesAesthetics(
                chevron_height_ratio=0.2,
                chevron_vertical_offset_ratio=1.0,
            ),
        )
        gene = pd.Series(
            {
                "start": 1000,
                "end": 6000,
                "strand": "+",
                "block_starts": [0, 4000],
                "block_sizes": [400, 500],
                "block_count": 2,
            }
        )

        ypos = 0.5
        genes._draw_gene_with_introns(ax, gene, ypos=ypos)

        first_chevron_y = _chevron_calls(ax)[0].args[1]
        chevron_peak = first_chevron_y[0]
        expected_peak = ypos - max(genes.interval_height * genes.chevron_height_ratio, 0.01)
        assert abs(chevron_peak - expected_peak) < 1e-6

    def test_draw_gene_with_introns_default_chevrons_symmetric_on_intron_line(self):
        ax = MagicMock()
        genes = Genes(data=pd.DataFrame())
        gene = pd.Series(
            {
                "start": 1000,
                "end": 6000,
                "strand": "+",
                "block_starts": [0, 4000],
                "block_sizes": [400, 500],
                "block_count": 2,
            }
        )

        ypos = 0.5
        genes._draw_gene_with_introns(ax, gene, ypos=ypos)

        first_chevron_y = _chevron_calls(ax)[0].args[1]

        assert abs(first_chevron_y[0] - ypos) < 1e-6

    def test_draw_gene_with_introns_reduces_chevron_width_when_gene_height_shrinks(self):
        fig, ax = plt.subplots(figsize=(6, 1.5))
        ax.set_xlim(1000, 6000)
        ax.set_ylim(0, 1)

        default_genes = Genes(data=pd.DataFrame(), aesthetics=GenesAesthetics(interval_height=0.2))
        short_genes = Genes(data=pd.DataFrame(), aesthetics=GenesAesthetics(interval_height=0.02))
        gene = pd.Series(
            {
                "start": 1000,
                "end": 6000,
                "strand": "+",
                "block_starts": [0, 4000],
                "block_sizes": [400, 500],
                "block_count": 2,
            }
        )

        default_genes._draw_gene_with_introns(ax, gene, ypos=0.7)
        default_marker_size = ax.lines[1].get_markersize()

        short_genes._draw_gene_with_introns(ax, gene, ypos=0.3)
        short_marker = ax.lines[3]
        short_marker_size = short_marker.get_markersize()

        plt.close(fig)

        assert short_marker.get_marker() == matplotlib.markers.CARETRIGHTBASE
        assert short_marker_size < default_marker_size
        assert short_marker_size >= 4.0

    def test_draw_gene_with_introns_keeps_chevrons_within_intron_bounds(self):
        fig, ax = plt.subplots(figsize=(6, 1.5))
        ax.set_xlim(1000, 6000)
        ax.set_ylim(0, 1)

        genes = Genes(data=pd.DataFrame())
        gene = pd.Series(
            {
                "start": 1000,
                "end": 6000,
                "strand": "+",
                "block_starts": [0, 4000],
                "block_sizes": [400, 500],
                "block_count": 2,
            }
        )

        intron_start = 1400
        intron_end = 5000
        genes._draw_gene_with_introns(ax, gene, ypos=0.5)

        chevrons = ax.lines[1:]
        assert chevrons

        for chevron in chevrons:
            assert chevron.get_clip_path() is not None
            left_extent, right_extent = genes._resolve_chevron_marker_extents(
                ax=ax,
                marker=chevron.get_marker(),
                marker_size=chevron.get_markersize(),
                marker_edge_width=chevron.get_markeredgewidth(),
            )
            boundary_padding = genes._resolve_chevron_boundary_padding(
                ax=ax,
                marker_edge_width=chevron.get_markeredgewidth(),
            )
            center = chevron.get_xdata()[0]
            assert center - left_extent - boundary_padding >= intron_start
            assert center + right_extent + boundary_padding <= intron_end

        plt.close(fig)

    def test_plot_genes_stagger_labels_uses_alternating_vertical_offsets(self, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1110],
                "end": [1120, 1130],
                "geneid": ["GENE_ALPHA", "GENE_BETA"],
                "block_count": [1, 1],
                "block_sizes": [[20], [20]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(display="collapsed", label_overlap_strategy="stagger"),
        )

        genes.plot_genes(mock_ax, genomic_region)

        assert len(mock_ax.text.call_args_list) == 2
        y_values = [call.args[1] for call in mock_ax.text.call_args_list]
        assert y_values[0] != y_values[1]

    def test_plot_genes_suppress_hides_overlapping_labels(self, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1110],
                "end": [1120, 1130],
                "geneid": ["GENE_ALPHA", "GENE_BETA"],
                "block_count": [1, 1],
                "block_sizes": [[20], [20]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(display="collapsed", label_overlap_strategy="suppress"),
        )

        genes.plot_genes(mock_ax, genomic_region)

        assert len(mock_ax.text.call_args_list) == 1

    def test_plot_genes_auto_strategy_shows_dense_labels(self, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1110],
                "end": [1120, 1130],
                "geneid": ["GENE_ALPHA", "GENE_BETA"],
                "block_count": [1, 1],
                "block_sizes": [[20], [20]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(data=data, aesthetics=GenesAesthetics(display="collapsed"))

        genes.plot_genes(mock_ax, genomic_region)

        assert len(mock_ax.text.call_args_list) == 2

    def test_plot_genes_default_enum_auto_strategy_stacks_dense_labels(
        self, mock_ax, genomic_region
    ):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1110],
                "end": [1120, 1130],
                "geneid": ["GENE_ALPHA", "GENE_BETA"],
                "block_count": [1, 1],
                "block_sizes": [[20], [20]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(data=data)

        genes.plot_genes(mock_ax, genomic_region)

        y_values = [call.args[1] for call in mock_ax.text.call_args_list]
        assert len(y_values) == 2
        assert len(set(y_values)) == 2

    def test_plot_genes_auto_strategy_avoids_connectors_by_default(
        self, mock_ax, genomic_region
    ):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1110],
                "end": [1120, 1130],
                "geneid": ["GENE_ALPHA", "GENE_BETA"],
                "block_count": [1, 1],
                "block_sizes": [[20], [20]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(display="collapsed", label_overlap_strategy="auto"),
        )

        genes.plot_genes(mock_ax, genomic_region)

        connector_calls = [
            call
            for call in mock_ax.plot.call_args_list
            if len(call.args) >= 2 and len(call.args[0]) == 2 and len(call.args[1]) == 2
        ]
        assert len(connector_calls) == 0

    def test_plot_genes_auto_strategy_can_draw_connectors_when_enabled(
        self, mock_ax, genomic_region
    ):
        data = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1940],
                "end": [1998],
                "geneid": ["GENE_AT_EDGE_WITH_LONG_NAME"],
                "block_count": [1],
                "block_sizes": [[58]],
                "block_starts": [[0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(
                display="collapsed",
                label_overlap_strategy="auto",
                label_connectors=True,
                label_max_chars=26,
            ),
        )

        genes.plot_genes(mock_ax, genomic_region)

        connector_calls = [
            call
            for call in mock_ax.plot.call_args_list
            if len(call.args) >= 2 and len(call.args[0]) == 2 and len(call.args[1]) == 2
        ]
        assert connector_calls

    def test_plot_genes_show_labels_false_skips_label_text(self, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1100],
                "end": [1200],
                "geneid": ["GENE_ALPHA"],
                "block_count": [1],
                "block_sizes": [[100]],
                "block_starts": [[0]],
            }
        )
        genes = Genes(data=data, aesthetics=GenesAesthetics(show_labels=False))

        genes.plot_genes(mock_ax, genomic_region)

        assert len(mock_ax.text.call_args_list) == 0

    def test_plot_genes_skips_labels_for_tiny_edge_overlap(self, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1990],
                "end": [3000],
                "geneid": ["EDGE_GENE"],
                "block_count": [1],
                "block_sizes": [[1010]],
                "block_starts": [[0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(
                display="collapsed",
                label_overlap_strategy="auto",
                label_min_overlap_bp=120,
                label_min_overlap_fraction=0.2,
            ),
        )

        genes.plot_genes(mock_ax, genomic_region)

        assert len(mock_ax.text.call_args_list) == 0

    def test_plot_genes_skips_labels_for_genes_outside_viewport(self, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [50, 1200],
                "end": [120, 1300],
                "geneid": ["OUTSIDE", "INSIDE"],
                "block_count": [1, 1],
                "block_sizes": [[70], [100]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(display="collapsed", label_overlap_strategy="auto"),
        )

        genes.plot_genes(mock_ax, genomic_region)

        labels = [call.args[2] for call in mock_ax.text.call_args_list]
        assert "OUTSIDE" not in labels
        assert "INSIDE" in labels

    def test_plot_genes_label_is_centered_on_gene_anchor(self, mock_ax, genomic_region):
        data = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1900],
                "end": [1990],
                "geneid": ["GENE_AT_RIGHT_EDGE"],
                "block_count": [1],
                "block_sizes": [[90]],
                "block_starts": [[0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(display="collapsed", label_overlap_strategy="suppress"),
        )

        genes.plot_genes(mock_ax, genomic_region)

        assert len(mock_ax.text.call_args_list) == 1
        label_x = mock_ax.text.call_args.args[0]
        expected_center = (1900 + 1990) / 2
        assert abs(label_x - expected_center) < 1e-6
        assert mock_ax.text.call_args.kwargs["ha"] == "center"

    def test_plot_genes_smart_expanded_auto_grows_rows_past_max(
        self, mock_ax, genomic_region
    ):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1110],
                "end": [1400, 1410],
                "geneid": ["GENE_ALPHA", "GENE_BETA"],
                "block_count": [1, 1],
                "block_sizes": [[300], [300]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(
                display="expanded",
                label_overlap_strategy="smart",
                max_number_of_rows=1,
            ),
        )

        genes.plot_genes(mock_ax, genomic_region)

        assert mock_ax.set_ylim.call_args.args[1] > 1.1

    def test_plot_genes_auto_expand_switches_from_collapsed_on_collisions(
        self, mock_ax, genomic_region
    ):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1110],
                "end": [1120, 1130],
                "geneid": ["GENE_ALPHA", "GENE_BETA"],
                "block_count": [1, 1],
                "block_sizes": [[20], [20]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(display="collapsed", label_overlap_strategy="auto_expand"),
        )

        genes.plot_genes(mock_ax, genomic_region)

        assert len(mock_ax.text.call_args_list) >= 1
        assert mock_ax.set_ylim.call_args.args[1] > 1.1

    def test_plot_genes_expanded_row_allocation_accounts_for_label_footprint(
        self, mock_ax, genomic_region
    ):
        data = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1200, 1260],
                "end": [1250, 1300],
                "geneid": ["VERY_LONG_GENE_LABEL_A", "VERY_LONG_GENE_LABEL_B"],
                "block_count": [1, 1],
                "block_sizes": [[50], [40]],
                "block_starts": [[0], [0]],
            }
        )
        genes = Genes(
            data=data,
            aesthetics=GenesAesthetics(display="expanded", label_overlap_strategy="suppress"),
        )

        genes.plot_genes(mock_ax, genomic_region)

        y_values = [call.args[1] for call in mock_ax.text.call_args_list]
        assert len(y_values) == 2
        assert len(set(y_values)) == 2
