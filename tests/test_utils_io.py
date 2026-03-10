import pandas as pd
from unittest.mock import MagicMock, patch

from plotnado.tracks.utils import (
    read_bed_regions,
    read_gtf_regions,
    intervals_to_dataframe,
)


class TestTrackIOUtils:
    def test_intervals_to_dataframe(self):
        records = [(10, 20, "name1"), (30, 40, "name2")]
        df = intervals_to_dataframe(records, "chr1")

        assert list(df["chrom"]) == ["chr1", "chr1"]
        assert list(df["start"]) == [10, 30]
        assert list(df["end"]) == [20, 40]

    @patch("plotnado.tracks.utils.pybigtools.open")
    def test_read_bed_regions_bigbed(self, mock_open):
        handle = MagicMock()
        handle.records.return_value = [(100, 200, "x"), (300, 400, "y")]
        mock_open.return_value.__enter__.return_value = handle

        df = read_bed_regions("mock.bb", "chr1", 50, 500)

        assert df.shape[0] == 2
        assert "chrom" in df.columns
        handle.records.assert_called_once_with("chr1", 50, 500)

    @patch("plotnado.tracks.utils._is_bigbed", return_value=False)
    @patch("plotnado.tracks.utils._read_bed_cached")
    def test_read_bed_regions_plain_bed(self, mock_cached, _mock_bigbed):
        mock_cached.return_value = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2"],
                "start": [100, 300, 100],
                "end": [200, 400, 200],
            }
        )

        df = read_bed_regions("mock.bed", "chr1", 150, 450)

        assert df.shape[0] == 2
        assert set(df["chrom"]) == {"chr1"}

    @patch("plotnado.tracks.utils._read_gtf_cached")
    def test_read_gtf_regions(self, mock_cached):
        mock_cached.return_value = pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [100, 300],
                "End": [200, 400],
                "Feature": ["exon", "exon"],
            }
        )

        df = read_gtf_regions("mock.gtf", "chr1", 50, 250)

        assert df.shape[0] == 1
        assert df.iloc[0]["start"] == 100
