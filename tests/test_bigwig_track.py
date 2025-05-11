"""
Tests for the BigWigTrack class.
"""

import pandas as pd
import pytest
import matplotlib.pyplot as plt
import numpy as np
from unittest.mock import patch, MagicMock

from plotnado.tracks import (
    BigWigTrack, 
    BigwigAesthetics, 
    GenomicRegion, 
    BedgraphDataFrame
)


class TestBigWigTrack:
    """Test cases for BigWigTrack class."""

    def test_init(self):
        """Test initializing a BigWigTrack instance."""
        aesthetics = BigwigAesthetics(color="red", fill=True)
        track = BigWigTrack(title="Test BigWig", aesthetics=aesthetics)
        
        assert track.title == "Test BigWig"
        assert track.aesthetics.color == "red"
        assert track.aesthetics.fill is True
        assert track.y_min is None
        assert track.y_max is None

    def test_fetch_from_df(self, genomic_region, bedgraph_df):
        """Test fetching data from a DataFrame."""
        # Setup
        track = BigWigTrack(
            title="Test BigWig",
            data=bedgraph_df,
            aesthetics=BigwigAesthetics()
        )
        
        # Test
        result = track._fetch_from_df(genomic_region)
        
        # Extract values in the genomic region
        expected = bedgraph_df[
            (bedgraph_df["chrom"] == genomic_region.chromosome) &
            (bedgraph_df["start"] >= genomic_region.start) &
            (bedgraph_df["end"] <= genomic_region.end)
        ]
        
        # Validate
        # If the region doesn't overlap with data, will return empty DataFrame
        if not expected.empty:
            assert result.shape[0] == expected.shape[0]

    def test_fetch_from_disk(self, genomic_region, test_bw_file):
        """Test fetching data from disk."""
        # Create track
        track = BigWigTrack(
            title="Test BigWig",
            data=test_bw_file,
            aesthetics=BigwigAesthetics()
        )
        
        # Test
        result = track._fetch_from_disk(genomic_region)

        # Extract values in the genomic region using real data  
        # This is a full test, so we need to ensure the data is correct
        # Validate the data
        assert result.shape[0] > 0
        assert all(result["chrom"] == genomic_region.chromosome)
        assert all(result["start"] >= genomic_region.start)
        assert all(result["end"] <= genomic_region.end)

    @patch.object(BigWigTrack, '_fetch_from_disk')
    @patch.object(BigWigTrack, '_fetch_from_df')
    def test_fetch_data_dataframe(self, mock_fetch_df, mock_fetch_disk, genomic_region, bedgraph_df):
        """Test fetch_data with DataFrame input."""
        # Setup
        mock_fetch_df.return_value = BedgraphDataFrame({
            'chrom': ['chr1'],
            'start': [100],
            'end': [200],
            'value': [0.5]
        })
        
        track = BigWigTrack(
            title="Test BigWig", 
            data=bedgraph_df,
            aesthetics=BigwigAesthetics()
        )
        
        # Test
        result = track.fetch_data(genomic_region)
        
        # Validate
        mock_fetch_disk.assert_not_called()
        mock_fetch_df.assert_called_once_with(genomic_region)

    def test_plot_stairs(self, mock_ax, genomic_region):
        """Test plotting stairs."""
        # Setup
        values = BedgraphDataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [100, 200],
            'end': [200, 300],
            'value': [0.5, 0.6]
        })
        
        aesthetics = BigwigAesthetics(
            color="red",
            fill=True,
            alpha=0.7,
            linewidth=2.0
        )
        
        track = BigWigTrack(
            title="Test BigWig",
            aesthetics=aesthetics
        )
        
        # Test
        track._plot_stairs(mock_ax, genomic_region, values)
        
        # Validate - check that the correct plotting method was called
        # Note: stairs returns an artist, so we check the number of calls
        assert mock_ax.stairs.call_count == 1
        args, kwargs = mock_ax.stairs.call_args
        assert kwargs['color'] == 'red'
        assert kwargs['alpha'] == 0.7
        assert kwargs['fill'] is True
        assert kwargs['linewidth'] == 2.0

    def test_plot_scatter(self, mock_ax, genomic_region):
        """Test plotting scatter."""
        # Setup
        values = BedgraphDataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [100, 200],
            'end': [200, 300],
            'value': [0.5, 0.6]
        })
        
        aesthetics = BigwigAesthetics(
            color="blue",
            alpha=0.7,
            scatter_point_size=5.0
        )
        
        track = BigWigTrack(
            title="Test BigWig",
            aesthetics=aesthetics
        )
        
        # Test
        track._plot_scatter(mock_ax, genomic_region, values)
        
        # Validate
        mock_ax.scatter.assert_called_once()
        args, kwargs = mock_ax.scatter.call_args
        assert kwargs['color'] == 'blue'
        assert kwargs['alpha'] == 0.7
        assert kwargs['s'] == 5.0

    @patch.object(BigWigTrack, 'fetch_data')
    @patch.object(BigWigTrack, '_plot_stairs')
    @patch.object(BigWigTrack, '_plot_scatter')
    def test_plot_std_style(self, mock_plot_scatter, mock_plot_stairs, mock_fetch_data, 
                         mock_ax, genomic_region):
        """Test plot with standard style."""
        # Setup
        mock_fetch_data.return_value = BedgraphDataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [100, 200],
            'end': [200, 300],
            'value': [0.5, 0.6]
        })
        
        aesthetics = BigwigAesthetics(
            style="std",
            min_value=0.0,
            max_value=1.0
        )
        
        track = BigWigTrack(
            title="Test BigWig",
            aesthetics=aesthetics
        )
        
        # Test
        track.plot(genomic_region, mock_ax)
        
        # Validate
        mock_fetch_data.assert_called_once_with(genomic_region)
        mock_plot_stairs.assert_called_once()
        mock_plot_scatter.assert_not_called()
        
        # Check that limits were set correctly
        mock_ax.set_xlim.assert_called_with(genomic_region.start, genomic_region.end)
        mock_ax.set_ylim.assert_called_with(ymin=0.0, ymax=1.0)
        
    @patch.object(BigWigTrack, 'fetch_data')
    @patch.object(BigWigTrack, '_plot_stairs')
    @patch.object(BigWigTrack, '_plot_scatter')
    def test_plot_scatter_style(self, mock_plot_scatter, mock_plot_stairs, mock_fetch_data, 
                             mock_ax, genomic_region):
        """Test plot with scatter style."""
        # Setup
        mock_fetch_data.return_value = BedgraphDataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [100, 200],
            'end': [200, 300],
            'value': [0.5, 0.6]
        })
        
        aesthetics = BigwigAesthetics(
            style="scatter",
            min_value=0.0,
            max_value=1.0
        )
        
        track = BigWigTrack(
            title="Test BigWig",
            aesthetics=aesthetics
        )
        
        # Test
        track.plot(genomic_region, mock_ax)
        
        # Validate
        mock_fetch_data.assert_called_once_with(genomic_region)
        mock_plot_scatter.assert_called_once()
        mock_plot_stairs.assert_not_called()
        
        # Check that limits were set correctly
        mock_ax.set_xlim.assert_called_with(genomic_region.start, genomic_region.end)
        mock_ax.set_ylim.assert_called_with(ymin=0.0, ymax=1.0)
