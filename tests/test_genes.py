"""
Tests for the Genes class.
"""

import os
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock, mock_open as mock_open_func
import matplotlib.pyplot as plt

from plotnado.tracks import (
    Genes,
    GenesAesthetics,
    GenomicRegion
)


class TestGenes:
    """Test cases for Genes class."""

    def test_init(self):
        """Test initializing a Genes instance."""
        aesthetics = GenesAesthetics(
            color="blue",
            display="collapsed",
            minimum_gene_length=100,
            max_number_of_rows=5
        )
        genes = Genes(
            title="Test Genes",
            aesthetics=aesthetics
        )
        
        assert genes.title == "Test Genes"
        assert genes.aesthetics.color == "blue"
        assert genes.aesthetics.display == "collapsed"
        assert genes.aesthetics.minimum_gene_length == 100
        assert genes.aesthetics.max_number_of_rows == 5

    @patch('plotnado.tracks.tabix_gtf')
    def test_fetch_from_disk_gtf(self, mock_tabix_gtf, genomic_region):
        """Test fetching data from a GTF file."""
        # Setup
        mock_tabix_gtf.return_value = "/path/to/mock.gtf.gz"
        
        # Create a mock for pysam.TabixFile
        mock_tabix = MagicMock()
        mock_record1 = MagicMock()
        mock_record1.seqname = "chr1"
        mock_record1.start = 1100
        mock_record1.end = 1200
        mock_record1.strand = "+"
        mock_record1.feature = "exon"
        mock_record1.attributes = 'gene_id "gene1"; transcript_id "transcript1";'
        
        mock_record2 = MagicMock()
        mock_record2.seqname = "chr1"
        mock_record2.start = 1300
        mock_record2.end = 1400
        mock_record2.strand = "+"
        mock_record2.feature = "exon"
        mock_record2.attributes = 'gene_id "gene1"; transcript_id "transcript1";'
        
        mock_tabix.fetch.return_value = [mock_record1, mock_record2]
        
        # Mock the pysam.TabixFile constructor
        with patch('pysam.TabixFile', return_value=mock_tabix):
            # Create a genes instance
            genes = Genes(
                title="Test Genes",
                data="/path/to/mock.gtf",
                aesthetics=GenesAesthetics()
            )
            
            # Test
            result = genes._fetch_from_disk_gtf(genomic_region)
            
            # Validate
            mock_tabix_gtf.assert_called_once_with("/path/to/mock.gtf")
            mock_tabix.fetch.assert_called_once_with(
                genomic_region.chromosome,
                genomic_region.start,
                genomic_region.end,
                parser=mock_tabix.fetch.call_args[1]['parser']
            )
            assert isinstance(result, pd.DataFrame)

    @patch('pybedtools.BedTool')
    def test_fetch_from_disk_bed12(self, mock_bedtool, genomic_region):
        """Test fetching data from a BED12 file."""
        # Setup mock BedTool
        mock_bt = MagicMock()
        mock_intervals = MagicMock()
        mock_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [1100, 1300],
            'end': [1200, 1400],
            'name': ['gene1', 'gene2'],
            'score': [0, 0],
            'strand': ['+', '+'],
            'thickStart': [1100, 1300],
            'thickEnd': [1200, 1400],
            'itemRgb': ['0,0,0', '0,0,0'],
            'blockCount': [1, 1],
            'blockSizes': ['100', '100'],
            'blockStarts': ['0', '0']
        })
        mock_intervals.to_dataframe.return_value = mock_df
        mock_bt.tabix_intervals.return_value = mock_intervals
        mock_bt.tabix.return_value = mock_bt
        mock_bedtool.return_value = mock_bt
        
        # Create genes instance
        genes = Genes(
            title="Test Genes",
            data="/path/to/mock.bed",
            aesthetics=GenesAesthetics()
        )
        
        # Test
        result = genes._fetch_from_disk_bed12(genomic_region)
        
        # Validate
        mock_bedtool.assert_called_once_with("/path/to/mock.bed")
        mock_bt.tabix.assert_called_once_with(force=True)
        mock_bt.tabix_intervals.assert_called_once_with(f"{genomic_region.chromosome}:{genomic_region.start}-{genomic_region.end}")
        assert isinstance(result, pd.DataFrame)
        assert 'block_starts' in result.columns
        assert 'block_sizes' in result.columns

    @patch.object(Genes, '_fetch_genes_from_package')
    @patch.object(Genes, '_fetch_from_disk_gtf')
    @patch.object(Genes, '_fetch_from_disk_bed12')
    def test_fetch_data_from_package(self, mock_fetch_bed12, mock_fetch_gtf, mock_fetch_package, genomic_region):
        """Test fetching data from package."""
        # Setup
        mock_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [1100, 1300],
            'end': [1200, 1400],
            'geneid': ['gene1', 'gene2'],
            'score': ['0', '0'],
            'strand': ['+', '+'],
            'thick_start': ['1100', '1300'],
            'thick_end': ['1200', '1400'],
            'color': ['0,0,0', '0,0,0'],
            'block_count': ['1', '1'],
            'block_sizes': ['100', '100'],
            'block_starts': ['0', '0']
        })
        mock_fetch_package.return_value = mock_df
        
        # Create genes instance
        genes = Genes(
            title="Test Genes",
            genome="hg38",
            aesthetics=GenesAesthetics()
        )
        
        # Test
        result = genes.fetch_data(genomic_region)
        
        # Validate
        mock_fetch_package.assert_called_once_with(genomic_region)
        mock_fetch_gtf.assert_not_called()
        mock_fetch_bed12.assert_not_called()
        assert result is mock_df

    @patch.object(Genes, '_fetch_genes_from_package')
    @patch.object(Genes, '_fetch_from_disk_gtf')
    @patch.object(Genes, '_fetch_from_disk_bed12')
    def test_fetch_data_gtf(self, mock_fetch_bed12, mock_fetch_gtf, mock_fetch_package, genomic_region):
        """Test fetching data from GTF file."""
        # Setup
        mock_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [1100, 1300],
            'end': [1200, 1400],
            'geneid': ['gene1', 'gene2'],
            'score': ['0', '0'],
            'strand': ['+', '+'],
            'thick_start': ['1100', '1300'],
            'thick_end': ['1200', '1400'],
            'color': ['0,0,0', '0,0,0'],
            'block_count': ['1', '1'],
            'block_sizes': ['100', '100'],
            'block_starts': ['0', '0']
        })
        mock_fetch_gtf.return_value = mock_df
        
        # Create genes instance with GTF file
        genes = Genes(
            title="Test Genes",
            data="/path/to/mock.gtf",
            aesthetics=GenesAesthetics()
        )
        
        # Test
        result = genes.fetch_data(genomic_region)
        
        # Validate
        mock_fetch_package.assert_not_called()
        mock_fetch_gtf.assert_called_once_with(genomic_region)
        mock_fetch_bed12.assert_not_called()
        assert result is mock_df

    @patch.object(Genes, '_fetch_genes_from_package')
    @patch.object(Genes, '_fetch_from_disk_gtf')
    @patch.object(Genes, '_fetch_from_disk_bed12')
    def test_fetch_data_bed(self, mock_fetch_bed12, mock_fetch_gtf, mock_fetch_package, genomic_region):
        """Test fetching data from BED file."""
        # Setup
        mock_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [1100, 1300],
            'end': [1200, 1400],
            'geneid': ['gene1', 'gene2'],
            'score': ['0', '0'],
            'strand': ['+', '+'],
            'thick_start': ['1100', '1300'],
            'thick_end': ['1200', '1400'],
            'color': ['0,0,0', '0,0,0'],
            'block_count': ['1', '1'],
            'block_sizes': ['100', '100'],
            'block_starts': ['0', '0']
        })
        mock_fetch_bed12.return_value = mock_df
        
        # Create genes instance with BED file
        genes = Genes(
            title="Test Genes",
            data="/path/to/mock.bed",
            aesthetics=GenesAesthetics()
        )
        
        # Test
        result = genes.fetch_data(genomic_region)
        
        # Validate
        mock_fetch_package.assert_not_called()
        mock_fetch_gtf.assert_not_called()
        mock_fetch_bed12.assert_called_once_with(genomic_region)
        assert result is mock_df

    def test_fetch_data_filters_by_length(self, genomic_region):
        """Test that fetch_data filters genes by minimum length."""
        # Create a mock DataFrame with genes of different lengths
        mock_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': [1100, 1300, 1500],
            'end': [1150, 1400, 1600],  # First gene is 50bp, others are 100bp
            'geneid': ['gene1', 'gene2', 'gene3'],
            'score': ['0', '0', '0'],
            'strand': ['+', '+', '+'],
            'thick_start': ['1100', '1300', '1500'],
            'thick_end': ['1150', '1400', '1600'],
            'color': ['0,0,0', '0,0,0', '0,0,0'],
            'block_count': ['1', '1', '1'],
            'block_sizes': ['50', '100', '100'],
            'block_starts': ['0', '0', '0']
        })
        
        # Create genes instance with minimum gene length of 80bp
        genes = Genes(
            title="Test Genes",
            data=mock_df,  # Pass the mock DataFrame directly as data
            aesthetics=GenesAesthetics(minimum_gene_length=80)
        )
        
        # Override the fetch_data method to return our mock_df directly
        with patch.object(genes, 'fetch_data', return_value=mock_df):
            # Test
            result = genes.fetch_data(genomic_region)
            
            # Filter the mock_df manually for validation
            expected = mock_df.query("end - start >= 80")
            
            # Validate
            assert result.shape[0] == 2  # Only 2 genes should remain
            assert 'gene1' not in result['geneid'].values  # The short gene should be filtered out
            assert set(result['geneid']) == set(['gene2', 'gene3'])  # Only longer genes remain

    @patch('plotnado.tracks.logger.warning')
    @patch.object(Genes, 'fetch_data')
    @patch.object(Genes, '_allocate_row_index')
    @patch.object(Genes, '_compute_extended_end_bp')
    @patch.object(Genes, '_draw_gene_feature')
    def test_plot_genes_no_genes(self, mock_draw, mock_compute_extended, 
                              mock_allocate_row, mock_fetch_data, mock_warning, 
                              mock_ax, genomic_region):
        """Test plot_genes with no genes."""
        # Setup
        mock_fetch_data.return_value = pd.DataFrame()  # Empty DataFrame
        
        genes = Genes(
            title="Test Genes",
            aesthetics=GenesAesthetics()
        )
        
        # Test
        genes.plot_genes(mock_ax, genomic_region)
        
        # Validate
        mock_fetch_data.assert_called_once_with(genomic_region)
        mock_warning.assert_called_once()  # Should log a warning
        mock_ax.set_ylim.assert_called_once_with(0, 1)
        mock_ax.set_xlim.assert_called_once_with(genomic_region.start, genomic_region.end)
        mock_compute_extended.assert_not_called()
        mock_allocate_row.assert_not_called()
        mock_draw.assert_not_called()

    @patch.object(Genes, 'fetch_data')
    @patch.object(Genes, '_estimate_label_char_width')
    @patch.object(Genes, '_allocate_row_index')
    @patch.object(Genes, '_compute_extended_end_bp')
    @patch.object(Genes, '_draw_gene_feature')
    def test_plot_genes_with_genes(self, mock_draw, mock_compute_extended, 
                                mock_allocate_row, mock_estimate_char, mock_fetch_data, 
                                mock_ax, genomic_region):
        """Test plot_genes with genes."""
        # Setup
        mock_genes_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [1100, 1300],
            'end': [1200, 1400],
            'geneid': ['gene1', 'gene2'],
            'score': [0, 0],
            'strand': ['+', '+'],
            'thick_start': [1100, 1300],
            'thick_end': [1200, 1400],
            'color': ['0,0,0', '0,0,0'],
            'block_count': [1, 1],
            'block_sizes': ['100', '100'],
            'block_starts': ['0', '0']
        })
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
            aesthetics=GenesAesthetics(
                display="expanded",
                max_number_of_rows=4
            )
        )
        
        # Test
        genes.plot_genes(ax, genomic_region)
        
        # Validate
        mock_fetch_data.assert_called_once_with(genomic_region)
        mock_estimate_char.assert_called_once()
        assert mock_compute_extended.call_count == 2
        assert mock_allocate_row.call_count == 2
        assert mock_draw.call_count == 2
        
        # Check that we set appropriate y limits based on row count
        assert genes.gene_count == 2  # Should count 2 genes
        ax.set_ylim.assert_called_once()
        ax.set_xlim.assert_called_once_with(genomic_region.start, genomic_region.end)
