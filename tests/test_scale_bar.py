"""
Tests for the ScaleBar class.
"""

import os
import pytest
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from unittest.mock import patch, call, MagicMock

from plotnado.tracks import (
    ScaleBar,
    ScaleBarAesthetics,
    GenomicRegion,
    get_human_readable_number_of_bp
)


class TestScaleBar:
    """Test cases for ScaleBar class."""

    def test_init(self):
        """Test initializing a ScaleBar instance."""
        aesthetics = ScaleBarAesthetics(
            color="blue",
            position="center",
            scale_distance=500
        )
        scale_bar = ScaleBar(
            title="Test Scale",
            aesthetics=aesthetics
        )
        
        assert scale_bar.title == "Test Scale"
        assert scale_bar.aesthetics.color == "blue"
        assert scale_bar.aesthetics.position == "center"
        assert scale_bar.aesthetics.scale_distance == 500

    def test_get_appropriate_scale(self):
        """Test the _get_appropriate_scale method."""
        scale_bar = ScaleBar(
            title="Test Scale",
            aesthetics=ScaleBarAesthetics()
        )
        
        # Test various lengths
        assert scale_bar._get_appropriate_scale(100) == 10
        assert scale_bar._get_appropriate_scale(1000) == 200
        assert scale_bar._get_appropriate_scale(10000) == 2000
        assert scale_bar._get_appropriate_scale(100000) == 20000
        
        # Test edge cases
        with pytest.raises(ValueError):
            scale_bar._get_appropriate_scale(0)
        with pytest.raises(ValueError):
            scale_bar._get_appropriate_scale(-1)

    def test_plot_left_position(self, mock_ax, genomic_region):
        """Test plotting ScaleBar with left position."""
        # Setup
        scale_bar = ScaleBar(
            title="Test Scale",
            aesthetics=ScaleBarAesthetics(
                position="left",
                scale_distance=200
            )
        )
        genomic_region = GenomicRegion(
            chromosome="chr1",
            start=1000,
            end=2000,
            strand="+"
        )
        
        # Test
        scale_bar.plot(mock_ax, genomic_region)
        
        # Validate
        # Should call plot for main scale bar and two ticks
        assert mock_ax.plot.call_count == 3
        mock_ax.text.assert_called_once()  # Scale text
        mock_ax.set_xlim.assert_called_once_with(1000, 2000)
        mock_ax.set_ylim.assert_called_once_with(0, 1)

    def test_plot_right_position(self, mock_ax, genomic_region):
        """Test plotting ScaleBar with right position."""
        # Setup
        scale_bar = ScaleBar(
            title="Test Scale",
            aesthetics=ScaleBarAesthetics(
                position="right",
                scale_distance=200
            )
        )
        
        # Test
        scale_bar.plot(mock_ax, genomic_region)
        
        # Validate
        # Should call plot for main scale bar and two ticks
        assert mock_ax.plot.call_count == 3
        mock_ax.text.assert_called_once()  # Scale text
        mock_ax.set_xlim.assert_called_once_with(1000, 2000)
        mock_ax.set_ylim.assert_called_once_with(0, 1)

    def test_plot_center_position(self, mock_ax, genomic_region):
        """Test plotting ScaleBar with center position."""
        # Setup
        scale_bar = ScaleBar(
            title="Test Scale",
            aesthetics=ScaleBarAesthetics(
                position="center",
                scale_distance=200
            )
        )
        
        # Test
        scale_bar.plot(mock_ax, genomic_region)
        
        # Validate
        # Should call plot for main scale bar and two ticks
        assert mock_ax.plot.call_count == 3
        mock_ax.text.assert_called_once()  # Scale text
        mock_ax.set_xlim.assert_called_once_with(1000, 2000)
        mock_ax.set_ylim.assert_called_once_with(0, 1)

    def test_plot_invalid_position(self, mock_ax, genomic_region):
        """Test plotting ScaleBar with invalid position."""
        # Setup
        scale_bar = ScaleBar(
            title="Test Scale",
            aesthetics=ScaleBarAesthetics(
                position="invalid",  # Invalid position
                scale_distance=200
            )
        )
        
        # Test & Validate
        with pytest.raises(ValueError, match="Position can only be"):
            scale_bar.plot(mock_ax, genomic_region)

    def test_plot_auto_scale_distance(self, mock_ax, genomic_region):
        """Test plotting ScaleBar with automatic scale distance."""
        # Setup
        scale_bar = ScaleBar(
            title="Test Scale",
            aesthetics=ScaleBarAesthetics(
                position="center",
                scale_distance=None  # Auto scale
            )
        )
        
        # Mock the _get_appropriate_scale method
        with patch.object(scale_bar, '_get_appropriate_scale', return_value=200) as mock_scale:
            # Test
            scale_bar.plot(mock_ax, genomic_region)
            
            # Validate
            mock_scale.assert_called_once_with(1000)  # 1000 is the length of genomic_region
            assert mock_ax.plot.call_count == 3
            mock_ax.text.assert_called_once()  # Scale text


class TestHumanReadableBp:
    """Test cases for get_human_readable_number_of_bp function."""
    
    def test_bp_scale(self):
        """Test bp scale."""
        assert get_human_readable_number_of_bp(100) == " 100 bp"
        assert get_human_readable_number_of_bp(500) == " 500 bp"
        assert get_human_readable_number_of_bp(999) == " 999 bp"
    
    def test_kb_scale(self):
        """Test kb scale."""
        assert get_human_readable_number_of_bp(1000) == " 1 kb"
        assert get_human_readable_number_of_bp(1500) == " 2 kb"
        assert get_human_readable_number_of_bp(10000) == " 10 kb"
        assert get_human_readable_number_of_bp(999999) == " 1000 kb"
    
    def test_mb_scale(self):
        """Test mb scale."""
        assert get_human_readable_number_of_bp(1000000) == " 1 mb"
        assert get_human_readable_number_of_bp(1500000) == " 2 mb"
        assert get_human_readable_number_of_bp(10000000) == " 10 mb"
