"""
Pytest configurations and fixtures for testing plotnado.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest
from matplotlib.testing.compare import compare_images

from plotnado.tracks import GenomicRegion


@pytest.fixture
def genomic_region():
    """Return a GenomicRegion instance for testing."""
    return GenomicRegion(chromosome="chr1", start=1000, end=2000, strand="+")


@pytest.fixture
def mock_ax():
    """Return a mocked matplotlib Axes instance for testing."""
    mock = MagicMock()

    # Add call_count attributes to plt methods for testing
    mock.text = MagicMock()
    mock.text.call_count = 0
    mock.text.assert_called_once = lambda: pytest.assume(mock.text.call_count == 1)
    mock.text.assert_called_once_with = (
        lambda *args, **kwargs: mock.text.assert_called_with(*args, **kwargs)
    )

    mock.plot = MagicMock()
    mock.plot.call_count = 0
    mock.plot.assert_called_once = lambda: pytest.assume(mock.plot.call_count == 1)
    mock.plot.assert_called_once_with = (
        lambda *args, **kwargs: mock.plot.assert_called_with(*args, **kwargs)
    )

    mock.set_xlim = MagicMock()
    mock.set_ylim = MagicMock()
    mock.set_xticks = MagicMock()
    mock.set_yticks = MagicMock()
    mock.spines = {
        "top": MagicMock(),
        "right": MagicMock(),
        "left": MagicMock(),
        "bottom": MagicMock(),
    }
    mock.xaxis = MagicMock()
    mock.yaxis = MagicMock()

    return mock


@pytest.fixture
def clean_up_temp_files():
    """Fixture to clean up temporary files after tests."""
    temp_files = []

    def _register_temp_file(filepath):
        temp_files.append(filepath)
        return filepath

    yield _register_temp_file

    # Clean up
    for filepath in temp_files:
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
            except Exception:
                pass


@pytest.fixture
def temp_bigwig(clean_up_temp_files):
    """Create a temporary bigwig file for testing."""
    import pyBigWig

    # Create a temporary bigwig file
    temp_bw = tempfile.NamedTemporaryFile(suffix=".bw", delete=False)
    temp_path = temp_bw.name
    temp_bw.close()

    # Write data to bigwig file
    bw = pyBigWig.open(temp_path, "w")
    bw.addHeader([("chr1", 10000)], maxZooms=0)

    # Add some test data
    start = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]
    end = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    bw.addEntries("chr1", start, ends=end, values=values)
    bw.close()

    # Register for cleanup
    clean_up_temp_files(temp_path)

    return temp_path


@pytest.fixture
def bedgraph_df():
    """Create a test bedgraph DataFrame."""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1", "chr1"],
            "start": [100, 200, 300, 400, 500],
            "end": [200, 300, 400, 500, 600],
            "value": [0.1, 0.2, 0.3, 0.4, 0.5],
        }
    )


@pytest.fixture
def bed_df():
    """Create a test BED DataFrame."""
    return pd.DataFrame(
        {
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [100, 300, 500],
            "End": [200, 400, 600],
            "Name": ["feature1", "feature2", "feature3"],
            "Score": [10, 20, 30],
            "Strand": ["+", "-", "+"],
        }
    )


@pytest.fixture
def temp_bed_file(bed_df, clean_up_temp_files):
    """Create a temporary BED file for testing."""
    temp_bed = tempfile.NamedTemporaryFile(suffix=".bed", delete=False)
    temp_path = temp_bed.name
    temp_bed.close()

    bed_df.to_csv(temp_path, sep="\t", index=False, header=False)

    # Register for cleanup
    clean_up_temp_files(temp_path)

    return temp_path


@pytest.fixture
def base_test_dir():
    """Return the base directory for test output."""
    base_dir = Path(os.path.dirname(os.path.abspath(__file__))) / "output"
    base_dir.mkdir(exist_ok=True)
    return base_dir


def assert_images_equal(actual_path, expected_path, tolerance=10):
    """Compare two images."""
    comparison = compare_images(actual_path, expected_path, tolerance)
    assert comparison is None, f"Images differ: {comparison}"
