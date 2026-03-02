import tempfile
import importlib.resources
import json
from pathlib import Path

import pandas as pd
import pytest

from plotnado import Figure
from plotnado.tracks import BigWigTrack, ScaleBar


class TestFigureRefactor:
    def test_extend_parameter(self):
        fig = Figure().add_track(ScaleBar())
        out = fig.plot("chr1:100-200", show=False, extend=0.5)

        xlim = out.axes[0].get_xlim()
        assert xlim[0] == 50
        assert xlim[1] == 250

    def test_plot_regions_list(self):
        fig = Figure().add_track(ScaleBar())
        plots = fig.plot_regions(["chr1:100-200", "chr1:300-400"], show=False)

        assert len(plots) == 2

    def test_plot_regions_grid_mode(self):
        fig = Figure().add_track(ScaleBar())
        plots = fig.plot_regions(
            ["chr1:100-200", "chr1:300-400", "chr1:500-600"],
            ncols=2,
            show=False,
        )

        assert len(plots) == 1
        assert len(plots[0].axes) == 4

    def test_plot_regions_invalid_ncols(self):
        fig = Figure().add_track(ScaleBar())

        with pytest.raises(ValueError, match="ncols must be >= 1"):
            fig.plot_regions(["chr1:100-200"], ncols=0, show=False)

    def test_plot_regions_bed_path(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as handle:
            handle.write("chr1\t100\t200\n")
            handle.write("chr1\t300\t400\n")
            bed_path = handle.name

        fig = Figure().add_track(ScaleBar())
        plots = fig.plot_regions(bed_path, show=False)

        assert len(plots) == 2
        Path(bed_path).unlink(missing_ok=True)

    def test_autoscale_group(self):
        df1 = pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [100, 150], "end": [150, 200], "value": [1.0, 2.0]}
        )
        df2 = pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [100, 150], "end": [150, 200], "value": [10.0, 20.0]}
        )

        fig = Figure()
        fig.add_track(BigWigTrack(data=df1, autoscale_group="g1"))
        fig.add_track(BigWigTrack(data=df2, autoscale_group="g1"))

        out = fig.plot("chr1:90-210", show=False)
        ylim1 = out.axes[0].get_ylim()
        ylim2 = out.axes[1].get_ylim()
        assert ylim1 == ylim2

    def test_toml_roundtrip(self):
        pytest.importorskip("tomli_w")
        fig = Figure().add_track(ScaleBar())

        with tempfile.NamedTemporaryFile(suffix=".toml", delete=False) as handle:
            path = handle.name

        fig.to_toml(path)
        loaded = Figure.from_toml(path)

        assert isinstance(loaded, Figure)
        assert len(loaded.tracks) == 1
        assert loaded.tracks[0].__class__.__name__ == "ScaleBar"
        Path(path).unlink(missing_ok=True)

    def test_toml_uses_track_tables(self):
        pytest.importorskip("tomli_w")
        fig = Figure().add_track(ScaleBar())

        with tempfile.NamedTemporaryFile(suffix=".toml", delete=False) as handle:
            path = handle.name

        fig.to_toml(path)
        content = Path(path).read_text()

        assert "[[tracks.ScaleBar]]" in content
        Path(path).unlink(missing_ok=True)

    def test_plot_gene_extend_math_with_bundled_hg38(self):
        bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
        with open(bed_prefix / "genes.json") as handle:
            mapping = json.load(handle)

        genes_df = pd.read_csv(bed_prefix / mapping["hg38"], sep="\t", header=None)
        genes_df.columns = [
            "chrom",
            "start",
            "end",
            "name",
            *[f"field_{i}" for i in range(max(0, genes_df.shape[1] - 4))],
        ]

        gene = str(genes_df.iloc[0]["name"])
        match = genes_df.loc[genes_df["name"].astype(str).str.upper() == gene.upper()].iloc[0]
        expected_length = int(match["end"]) - int(match["start"])
        expected_extend = int(expected_length * 0.5)

        fig = Figure().add_track(ScaleBar())
        out = fig.plot_gene(gene, extend=0.5, show=False)

        xlim = out.axes[0].get_xlim()
        assert xlim[0] == int(match["start"]) - expected_extend
        assert xlim[1] == int(match["end"]) + expected_extend

    def test_plot_gene_adds_genes_track_when_missing(self):
        fig = Figure().add_track(ScaleBar())
        out = fig.plot_gene("DDX11L1", show=False)

        assert out is not None
        assert len(out.axes) >= 2
