import pandas as pd

from plotnado.tracks import BedTrack, GenomicRegion


class _FakePyRanges:
    def __init__(self, df: pd.DataFrame):
        self.df = df


class TestBedTrackMemory:
    def test_dataframe_data(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2"],
                "start": [100, 300, 100],
                "end": [200, 400, 200],
            }
        )
        track = BedTrack(data=df)
        gr = GenomicRegion(chromosome="chr1", start=50, end=350)

        out = track.fetch_data(gr)
        assert out.shape[0] == 2

    def test_pyranges_like_data(self):
        df = pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [100, 300],
                "End": [200, 400],
            }
        )
        track = BedTrack(data=_FakePyRanges(df))
        gr = GenomicRegion(chromosome="chr1", start=50, end=350)

        out = track.fetch_data(gr)
        assert set(out.columns) >= {"chrom", "start", "end"}
        assert out.shape[0] == 2
