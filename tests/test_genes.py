"""
Tests for the Genes class.
"""

import os
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock, mock_open as mock_open_func
import matplotlib.pyplot as plt

from plotnado.tracks import Genes, GenesAesthetics, GenomicRegion


class TestGenes:
    """Test cases for Genes class."""

    def test_init(self):
        """Test initializing a Genes instance."""
        aesthetics = GenesAesthetics(
            color="blue",
            display="collapsed",
            minimum_gene_length=100,
            max_number_of_rows=5,
        )
        genes = Genes(title="Test Genes", aesthetics=aesthetics)

        assert genes.title == "Test Genes"
        assert genes.aesthetics.color == "blue"
        assert genes.aesthetics.display == "collapsed"
        assert genes.aesthetics.minimum_gene_length == 100
        assert genes.aesthetics.max_number_of_rows == 5

    def test_fetch_from_disk_gtf(self, test_gtf_file, genomic_region_chr21):
        """Test fetching data from a GTF file."""

        genes = Genes(aesthetics=GenesAesthetics(), data=test_gtf_file)

        genes_df = genes._fetch_from_disk_gtf(genomic_region_chr21)
        assert isinstance(genes_df, pd.DataFrame)
        assert not genes_df.empty
        assert genes_df.shape[0] == 1

    def test_fetch_from_disk_bed12(self, genomic_region_chr21, test_bed12_file):
        """Test fetching data from a BED12 file."""
        genes = Genes(aesthetics=GenesAesthetics(), data=test_bed12_file)
        genes_df = genes._fetch_from_disk_bed12(genomic_region_chr21)
        assert isinstance(genes_df, pd.DataFrame)
        assert not genes_df.empty


    @patch('os.access')
    @patch('os.path.dirname')
    def test_handle_write_permission_error(self, mock_dirname, mock_gettempdir, mock_access, genomic_region, test_gtf_file, mock_ax):
        """Test that Genes handles write permission errors."""
        # Setup
        mock_access.return_value = False
        mock_dirname.return_value = "/fake/readonly/dir"

        # Create a Genes instance with the test GTF file
        genes = Genes(
            title="Test Genes",
            aesthetics=GenesAesthetics(),
            data=test_gtf_file
        )

        # Check data fetching works despite persmission issues
        data = genes.fetch_data(genomic_region)
        assert not data.empty
        
        # Test the plotting functionality to ensure it works with write permission issues
        genes.plot_genes(ax=mock_ax, gr=genomic_region)
        
        # Verify expected behavior
        mock_access.assert_called()  # Permission check was performed
        mock_gettempdir.assert_called()  # Temp directory was used as fallback
        
        # Verify the plot was set up correctly
        mock_ax.set_xlim.assert_called_once_with(genomic_region.start, genomic_region.end)
        mock_ax.set_ylim.assert_called_once_with(0, 1)

            

    @patch.object(Genes, "fetch_data")
    @patch.object(Genes, "_allocate_row_index")
    @patch.object(Genes, "_compute_extended_end_bp")
    @patch.object(Genes, "_draw_gene_feature")
    def test_plot_genes_no_genes(
        self,
        mock_draw,
        mock_compute_extended,
        mock_allocate_row,
        mock_fetch_data,
        mock_warning,
        mock_ax,
        genomic_region,
    ):
        """Test plot_genes with no genes."""
        # Setup
        mock_fetch_data.return_value = pd.DataFrame()  # Empty DataFrame

        genes = Genes(title="Test Genes", aesthetics=GenesAesthetics())

        # Test
        genes.plot_genes(ax=mock_ax, gr=genomic_region)

        # Validate
        mock_fetch_data.assert_called_once_with(genomic_region)
        mock_ax.set_ylim.assert_called_once_with(0, 1)
        mock_ax.set_xlim.assert_called_once_with(
            genomic_region.start, genomic_region.end
        )
        mock_compute_extended.assert_not_called()
        mock_allocate_row.assert_not_called()
        mock_draw.assert_not_called()

    @patch.object(Genes, "fetch_data")
    @patch.object(Genes, "_estimate_label_char_width")
    @patch.object(Genes, "_allocate_row_index")
    @patch.object(Genes, "_compute_extended_end_bp")
    @patch.object(Genes, "_draw_gene_feature")
    def test_plot_genes_with_genes(
        self,
        mock_draw,
        mock_compute_extended,
        mock_allocate_row,
        mock_estimate_char,
        mock_fetch_data,
        mock_ax,
        genomic_region,
    ):
        """Test plot_genes with genes."""
        # Setup
        mock_genes_df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1100, 1300],
                "end": [1200, 1400],
                "geneid": ["gene1", "gene2"],
                "score": [0, 0],
                "strand": ["+", "+"],
                "thick_start": [1100, 1300],
                "thick_end": [1200, 1400],
                "color": ["0,0,0", "0,0,0"],
                "block_count": [1, 1],
                "block_sizes": ["100", "100"],
                "block_starts": ["0", "0"],
            }
        )
        mock_fetch_data.return_value = mock_genes_df

        # Mock the _compute_extended_end_bp method
        mock_compute_extended.side_effect = [1220, 1420]  # Add 20bp for label padding

        # Mock the _allocate_row_index method
        mock_allocate_row.side_effect = [0, 1]  # First gene row 0, second gene row 1

        # Create a figure to get a real figure width
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 5)

        # Create genes instance
        genes = Genes(
            title="Test Genes",
            aesthetics=GenesAesthetics(display="expanded", max_number_of_rows=4),
        )

        # Test
        genes.plot_genes(ax=mock_ax, gr=genomic_region)

        # Validate
        mock_fetch_data.assert_called_once_with(genomic_region)
        mock_estimate_char.assert_called_once()
        assert mock_compute_extended.call_count == 2
        assert mock_allocate_row.call_count == 2
        assert mock_draw.call_count == 2

        # Check that we set appropriate y limits based on row count
        assert genes.gene_count == 2  # Should count 2 genes
        mock_ax.set_ylim.assert_called_once()
        mock_ax.set_xlim.assert_called_once_with(
            genomic_region.start, genomic_region.end
        )
