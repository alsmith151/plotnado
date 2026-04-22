"""
Tests for the BedgraphDataFrame class.
"""

import pytest
import pandas as pd
from pandera.errors import SchemaError, SchemaErrors

from plotnado.tracks import BedgraphDataFrame


class TestBedgraphDataFrame:
    """Test cases for BedgraphDataFrame class."""

    def test_valid_creation(self):
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': [100, 200, 300],
            'end': [200, 300, 400],
            'value': [0.1, 0.2, 0.3],
        }
        df = BedgraphDataFrame(data)
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ['chrom', 'start', 'end', 'value']
        assert df.shape == (3, 4)

    def test_missing_column_raises(self):
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': [100, 200, 300],
            'end': [200, 300, 400],
        }
        with pytest.raises((SchemaError, SchemaErrors)):
            BedgraphDataFrame(data)

    def test_extra_column_raises(self):
        """Schema is strict=True so extra columns are rejected."""
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': [100, 200, 300],
            'end': [200, 300, 400],
            'value': [0.1, 0.2, 0.3],
            'extra': ['a', 'b', 'c'],
        }
        with pytest.raises(SchemaErrors):
            BedgraphDataFrame(data)

    def test_string_columns_are_coerced(self):
        data = {
            'chrom': ['chr1', 'chr1', 'chr1'],
            'start': ['100', '200', '300'],
            'end': ['200', '300', '400'],
            'value': ['0.1', '0.2', '0.3'],
        }
        df = BedgraphDataFrame(data)
        assert df['start'].dtype == 'int64'
        assert df['end'].dtype == 'int64'
        assert df['value'].dtype == 'float64'

    def test_unconvertible_value_raises(self):
        data = {
            'chrom': ['chr1'],
            'start': [100],
            'end': [200],
            'value': ['not_a_number'],
        }
        with pytest.raises((SchemaErrors, Exception)):
            BedgraphDataFrame(data)
