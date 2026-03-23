"""CLI integration tests via typer.testing.CliRunner."""

import pytest
from pathlib import Path
from typer.testing import CliRunner

from plotnado.cli.cli import app

runner = CliRunner(mix_stderr=False)


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
