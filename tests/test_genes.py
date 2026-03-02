"""Tests for the refactored Genes track."""

from unittest.mock import patch

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
