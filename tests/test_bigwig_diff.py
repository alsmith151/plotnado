import numpy as np
import pandas as pd
from unittest.mock import patch

from plotnado.tracks import BigWigDiff, GenomicRegion


class TestBigWigDiff:
    @staticmethod
    def _df_a() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1000, 1100],
                "end": [1100, 1200],
                "value": [4.0, 8.0],
            }
        )

    @staticmethod
    def _df_b() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1000, 1100],
                "end": [1100, 1200],
                "value": [2.0, 4.0],
            }
        )

    @patch("plotnado.tracks.bigwig.BigWigTrack.fetch_data")
    def test_fetch_subtract(self, mock_fetch):
        mock_fetch.side_effect = [self._df_a(), self._df_b()]
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)

        track = BigWigDiff(file_a="a.bw", file_b="b.bw", method="subtract")
        data = track.fetch_data(gr)

        assert np.allclose(data["value"].to_numpy(), np.array([2.0, 4.0]))

    @patch("plotnado.tracks.bigwig.BigWigTrack.fetch_data")
    def test_fetch_ratio(self, mock_fetch):
        mock_fetch.side_effect = [self._df_a(), self._df_b()]
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)

        track = BigWigDiff(file_a="a.bw", file_b="b.bw", method="ratio")
        data = track.fetch_data(gr)

        assert np.allclose(data["value"].to_numpy(), np.array([2.0, 2.0]))

    @patch.object(BigWigDiff, "fetch_data")
    def test_plot_empty(self, mock_fetch, mock_ax):
        mock_fetch.return_value = pd.DataFrame(columns=["x", "value"])
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)

        track = BigWigDiff(file_a="a.bw", file_b="b.bw")
        track.plot(mock_ax, gr)

        mock_ax.set_xlim.assert_called_once_with(gr.start, gr.end)
        mock_ax.set_ylim.assert_called_once_with(-1, 1)

    @patch.object(BigWigDiff, "fetch_data")
    def test_plot_uses_interval_bars(self, mock_fetch, mock_ax):
        mock_fetch.return_value = pd.DataFrame(
            {
                "x": [1050, 1150],
                "start": [1000, 1100],
                "end": [1100, 1200],
                "value": [2.0, -1.0],
            }
        )
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)

        track = BigWigDiff(file_a="a.bw", file_b="b.bw")
        track.plot(mock_ax, gr)

        assert mock_ax.bar.call_count == 2
        mock_ax.axhline.assert_called_once()
