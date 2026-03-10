import pandas as pd
from unittest.mock import patch

from plotnado.tracks import BigWigCollection, GenomicRegion


class TestBigWigCollection:
    @staticmethod
    def _df() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [1000, 1100],
                "end": [1100, 1200],
                "value": [1.0, 2.0],
            }
        )

    @patch("plotnado.tracks.bigwig.BigWigTrack.fetch_data")
    def test_fetch_data(self, mock_fetch):
        mock_fetch.return_value = self._df()
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)

        track = BigWigCollection(files=["a.bw", "b.bw"])
        result = track.fetch_data(gr)

        assert len(result) == 2
        assert mock_fetch.call_count == 2

    @patch("plotnado.tracks.bigwig.BigWigTrack.fetch_data")
    def test_plot_overlay(self, mock_fetch, mock_ax):
        mock_fetch.return_value = self._df()
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)

        track = BigWigCollection(files=["a.bw", "b.bw"])
        track.plot(mock_ax, gr)

        mock_ax.set_xlim.assert_called_with(gr.start, gr.end)
        assert mock_ax.set_ylim.call_count >= 1

    @patch("plotnado.tracks.bigwig.BigWigTrack.fetch_data")
    def test_plot_stacked(self, mock_fetch, mock_ax):
        mock_fetch.return_value = self._df()
        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)

        track = BigWigCollection(files=["a.bw", "b.bw"], aesthetics={"style": "stacked"})
        track.plot(mock_ax, gr)

        assert mock_ax.plot.call_count >= 1
