"""CLI integration tests via typer.testing.CliRunner."""

import inspect
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from typer.testing import CliRunner

from plotnado.cli.cli import app

_cli_runner_kwargs = {}
if "mix_stderr" in inspect.signature(CliRunner).parameters:
    _cli_runner_kwargs["mix_stderr"] = False
runner = CliRunner(**_cli_runner_kwargs)


# ---------------------------------------------------------------------------
# Fixture YAMLs
# ---------------------------------------------------------------------------

MINIMAL_YAML = """\
genome: hg38
tracks:
  - path: nonexistent.bw
    type: bigwig
    title: Track A
    group: g1
"""

YAML_WITH_BAD_GROUP_REF = """\
tracks:
  - path: nonexistent.bw
    type: bigwig
    title: Track A
groups:
  - name: g1
    tracks: ["WRONG TITLE"]
    autocolor: false
"""


@pytest.fixture
def minimal_yaml(tmp_path) -> Path:
    p = tmp_path / "template.yaml"
    p.write_text(MINIMAL_YAML)
    return p


@pytest.fixture
def bad_group_yaml(tmp_path) -> Path:
    p = tmp_path / "bad.yaml"
    p.write_text(YAML_WITH_BAD_GROUP_REF)
    return p


@pytest.fixture
def plot_template(tmp_path, test_bed12_file) -> tuple[Path, Path]:
    p = tmp_path / "plot.yaml"
    p.write_text(
        f"""\
genome: hg38
width: 14
track_height: 1.5
guides:
  genes: true
tracks:
  - path: {test_bed12_file}
    type: bed
    title: Signal
    group: sample1
    color: tomato
    height: 2.0
    options:
      data_range_style: text
groups:
  - name: sample1
    tracks: ["Signal"]
    autocolor: false
"""
    )
    return p, test_bed12_file


@pytest.fixture
def cooler_plot_template(tmp_path, test_cooler_file) -> tuple[Path, Path]:
    p = tmp_path / "cooler.yaml"
    p.write_text(
        f"""\
genome: hg38
width: 10
guides:
  genes: false
tracks:
  - path: {test_cooler_file}
    type: cooler
    title: Contacts
    height: 2.5
    options:
      cmap: magma
"""
    )
    return p, test_cooler_file


# ---------------------------------------------------------------------------
# init --auto
# ---------------------------------------------------------------------------

class TestInitAuto:
    def test_auto_with_bw_files(self, tmp_path):
        bw = tmp_path / "sample1_H3K27ac.bw"
        bw.touch()
        out = tmp_path / "out.yaml"
        result = runner.invoke(app, ["init", str(bw), "--auto", "--output", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_auto_generates_annotated_yaml(self, tmp_path):
        bw = tmp_path / "sample1_H3K27ac.bw"
        bw.touch()
        out = tmp_path / "out.yaml"
        runner.invoke(app, ["init", str(bw), "--auto", "--output", str(out)])
        content = out.read_text()
        assert content.startswith("# PlotNado template")

    def test_auto_track_ordering(self, tmp_path):
        """BigWig tracks should appear before narrowpeak in output."""
        bw = tmp_path / "sample1_H3K27ac.bw"
        np_ = tmp_path / "peaks.narrowpeak"
        bw.touch()
        np_.touch()
        out = tmp_path / "out.yaml"
        runner.invoke(app, ["init", str(np_), str(bw), "--auto", "--output", str(out)])

        import yaml
        data = yaml.safe_load(out.read_text())
        types = [t["type"] for t in data["tracks"]]
        bw_idx = types.index("bigwig")
        np_idx = types.index("narrowpeak")
        assert bw_idx < np_idx

    def test_color_assigned_in_output(self, tmp_path):
        bw = tmp_path / "sample1_H3K27ac.bw"
        bw.touch()
        out = tmp_path / "out.yaml"
        runner.invoke(app, ["init", str(bw), "--auto", "--output", str(out)])

        import yaml
        data = yaml.safe_load(out.read_text())
        assert data["tracks"][0].get("color") is not None


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

class TestValidate:
    def test_validates_minimal_yaml(self, minimal_yaml):
        result = runner.invoke(app, ["validate", str(minimal_yaml)])
        # File doesn't exist on disk → warns but still exits (may be 1 due to missing file)
        assert "loaded" in result.output.lower() or "Template" in result.output

    def test_catches_bad_group_ref(self, bad_group_yaml):
        result = runner.invoke(app, ["validate", str(bad_group_yaml)])
        assert result.exit_code == 1
        assert "Group reference error" in result.output or "not found" in result.output

    def test_file_not_found_error(self, tmp_path):
        result = runner.invoke(app, ["validate", str(tmp_path / "missing.yaml")])
        assert result.exit_code == 1
        assert "not found" in result.output.lower()

    def test_explain_flag_removed(self):
        """The --explain flag should no longer exist."""
        result = runner.invoke(app, ["validate", "--help"])
        assert "--explain" not in result.output


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------

class TestPlot:
    def test_plot_builds_figure_from_template_and_saves(self, plot_template, monkeypatch, tmp_path):
        plot_yaml, test_bed12_file = plot_template
        monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mplconfig"))

        out = tmp_path / "custom.png"
        result = runner.invoke(
            app,
            [
                "plot",
                str(plot_yaml),
                "--region",
                "chr1:1-10",
                "--output",
                str(out),
                "--width",
                "18",
                "--dpi",
                "300",
            ],
        )

        assert result.exit_code == 0, result.output
        assert out.exists()
        assert out.stat().st_size > 0
        assert "Compiled render plan" in result.output
        assert "Saved plot" in result.output

        image = plt.imread(out)
        assert image.size > 0
        assert image.shape[1] > image.shape[0]

    def test_plot_generates_one_output_per_region(self, plot_template, monkeypatch):
        plot_yaml, _ = plot_template
        monkeypatch.setenv("MPLCONFIGDIR", str(plot_yaml.parent / "mplconfig"))
        monkeypatch.chdir(plot_yaml.parent)

        result = runner.invoke(
            app,
            [
                "plot",
                str(plot_yaml),
                "--region",
                "chr21:5022531-5046683",
                "--region",
                "chr21:5030000-5040000",
            ],
        )

        assert result.exit_code == 0, result.output
        out1 = plot_yaml.parent / "plot_chr21_5022531_5046683_+_.png"
        out2 = plot_yaml.parent / "plot_chr21_5030000_5040000_+_.png"
        assert out1.exists()
        assert out2.exists()
        assert out1.stat().st_size > 0
        assert out2.stat().st_size > 0

    def test_plot_renders_cooler_template(self, cooler_plot_template, monkeypatch, tmp_path):
        plot_yaml, _ = cooler_plot_template
        monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mplconfig"))

        out = tmp_path / "cooler.png"
        result = runner.invoke(
            app,
            [
                "plot",
                str(plot_yaml),
                "--region",
                "chr1:2000000-10000000",
                "--output",
                str(out),
                "--dpi",
                "150",
            ],
        )

        assert result.exit_code == 0, result.output
        assert out.exists()
        assert out.stat().st_size > 0

        image = plt.imread(out)
        assert image.size > 0
        assert image.ndim == 3

    def test_plot_resolves_gene_names_from_template_genome(self, plot_template, monkeypatch, tmp_path):
        plot_yaml, _ = plot_template
        monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mplconfig"))

        out = tmp_path / "gnaq.png"
        result = runner.invoke(
            app,
            [
                "plot",
                str(plot_yaml),
                "--region",
                "GNAQ",
                "--output",
                str(out),
            ],
        )

        assert result.exit_code == 0, result.output
        assert "Resolved gene:" in result.output
        assert "chr9:77716097-78031811(+)" in result.output
        assert out.exists()
        assert out.stat().st_size > 0
