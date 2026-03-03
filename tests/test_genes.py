"""Tests for the refactored Genes track."""

from unittest.mock import MagicMock, patch

import pandas as pd

from plotnado.tracks import Genes, GenesAesthetics, GenomicRegion


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

        assert ax.plot.call_count >= 3

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

        first_chevron_x = ax.plot.call_args_list[1].args[0]
        assert first_chevron_x[0] > first_chevron_x[1]

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
        first_chevron_x = ax.plot.call_args_list[1].args[0]
        chevron_mid = sum(first_chevron_x) / 2

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

        chevron_line_calls = ax.plot.call_args_list[1::2]
        centers = [sum(call.args[0]) / 2 for call in chevron_line_calls]

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
        chevron_line_calls = ax.plot.call_args_list[1::2]
        centers = [sum(call.args[0]) / 2 for call in chevron_line_calls]

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

        chevron_line_calls = ax.plot.call_args_list[1::2]
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

        first_chevron_y = ax.plot.call_args_list[1].args[1]
        chevron_peak = first_chevron_y[1]
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

        lower_arm = ax.plot.call_args_list[1].args[1]
        upper_arm = ax.plot.call_args_list[2].args[1]

        assert abs(lower_arm[1] - ypos) < 1e-6
        assert abs(upper_arm[1] - ypos) < 1e-6
        assert abs((ypos - lower_arm[0]) - (upper_arm[0] - ypos)) < 1e-6
