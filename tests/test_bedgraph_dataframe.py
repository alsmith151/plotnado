"""
Tests for the BedgraphDataFrame class.
"""

import pytest
import pandas as pd
import pandera

from plotnado.tracks import BedgraphDataFrame


class TestBedgraphDataFrame:
    """Test cases for BedgraphDataFrame class."""

    def test_valid_creation(self):
        """Test creating a valid BedgraphDataFrame."""
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': [100, 200, 300],
            'end': [200, 300, 400],
            'value': [0.1, 0.2, 0.3]
        }
        # Use the function directly for testing
        df = BedgraphDataFrame(data)

        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ['chrom', 'start', 'end', 'value']
        assert df.shape == (3, 4)

    def test_missing_column(self):
        """Test creating a BedgraphDataFrame with a missing column."""
        # Missing 'value' column
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': [100, 200, 300],
            'end': [200, 300, 400]
        }
        
        with pytest.raises(pandera.errors.SchemaError):
            BedgraphDataFrame(data)

    def test_wrong_type(self):
        """Test creating a BedgraphDataFrame with wrong column types."""
        # 'start' column has strings instead of integers
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': ['100', '200', '300'],  # Strings instead of integers
            'end': [200, 300, 400],
            'value': [0.1, 0.2, 0.3]
        }
        
        # Should coerce strings to integers
        df = BedgraphDataFrame(data)
        assert df['start'].dtype == 'int64'

    def test_coercion(self):
        """Test that the coercion setting works."""
        # All columns have string values
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': ['100', '200', '300'],
            'end': ['200', '300', '400'],
            'value': ['0.1', '0.2', '0.3']
        }
        
        # Should coerce to appropriate types
        df = BedgraphDataFrame(data)
        assert df['chrom'].dtype == 'object'  # Strings
        assert df['start'].dtype == 'int64'   # Integers
        assert df['end'].dtype == 'int64'     # Integers
        assert df['value'].dtype == 'float64' # Floats
        
    def test_extra_columns(self):
        """Test creating a BedgraphDataFrame with extra columns."""
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': [100, 200, 300],
            'end': [200, 300, 400],
            'value': [0.1, 0.2, 0.3],
            'extra': ['a', 'b', 'c']  # Extra column
        }
        
        # Should raise error due to extra columns
        with pytest.raises(pandera.errors.SchemaError):
            BedgraphDataFrame(data)
