import sys
import types
import numpy as np

from plotnado.tracks import CoolerTrack, CoolerAverage, GenomicRegion


class _FakeMatrix:
    def __init__(self, matrix: np.ndarray):
        self._matrix = matrix

    def fetch(self, region: str):
        return self._matrix


class _FakeCooler:
    def __init__(self, uri: str):
        self.uri = uri

    def matrix(self, balance: bool = True):
        return _FakeMatrix(np.array([[1.0, 2.0], [3.0, 4.0]]))


class TestCoolerTrack:
    def test_fetch_data_log2(self, monkeypatch):
        fake_module = types.SimpleNamespace(Cooler=_FakeCooler)
        monkeypatch.setitem(sys.modules, "cooler", fake_module)

        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)
        track = CoolerTrack(file="a.cool", transform="log2")

        matrix = track.fetch_data(gr)
        assert matrix.shape == (2, 2)
        assert np.isclose(matrix[0, 0], 0.0)

    def test_uri_for_mcool_resolution(self):
        track = CoolerTrack(file="a.mcool", resolution=10000)
        assert track._cooler_uri() == "a.mcool::resolutions/10000"

    def test_average_track(self, monkeypatch):
        fake_module = types.SimpleNamespace(Cooler=_FakeCooler)
        monkeypatch.setitem(sys.modules, "cooler", fake_module)

        gr = GenomicRegion(chromosome="chr1", start=1000, end=2000)
        track = CoolerAverage(files=["a.cool", "b.cool"])

        matrix = track.fetch_data(gr)
        assert matrix.shape == (2, 2)
        assert np.allclose(matrix, np.array([[1.0, 2.0], [3.0, 4.0]]))
