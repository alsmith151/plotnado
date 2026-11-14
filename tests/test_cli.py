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


@pytest.mark.parametrize("args", [["--version"], ["-V"], ["version"]])
def test_version_reports_installed_version(args):
    from plotnado.cli.cli import get_version

    result = runner.invoke(app, args)
    assert result.exit_code == 0
    assert result.output.strip() == f"plotnado {get_version()}"


def _run_cli(*args, cwd=None, input=None):
    """Run the real CLI in a subprocess (needed where Python's own stdio/exit code matter)."""
    import subprocess
    import sys

    return subprocess.run(
        [sys.executable, "-c", "from plotnado.cli.cli import main; main()", *args],
        capture_output=True,
        text=True,
        cwd=cwd,
        input=input,
    )


def test_help_lists_environment_commands():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    for command in ("version", "python", "kernel"):
        assert command in result.output


def test_python_command_runs_script_with_args_and_exit_code(tmp_path):
    script = tmp_path / "s.py"
    script.write_text("import sys, plotnado\nprint('ok', sys.argv[1:])\nsys.exit(3)\n")
    result = _run_cli("python", str(script), "-x", "--flag")
    assert result.returncode == 3, result.stderr
    assert "ok ['-x', '--flag']" in result.stdout


def test_python_command_runs_inline_code():
    result = _run_cli("python", "-c", "import plotnado; print('version', plotnado.__version__)")
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("version ")


def test_python_command_runs_module(tmp_path):
    pkg = tmp_path / "mypkg"
    pkg.mkdir()
    (pkg / "__init__.py").write_text("")
    (pkg / "plots.py").write_text("import plotnado\nprint('module ran')\n")
    result = _run_cli("python", "-m", "mypkg.plots", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert "module ran" in result.stdout


def test_python_command_reads_stdin_without_args():
    result = _run_cli("python", input="import plotnado\nprint('from stdin')\n")
    assert result.returncode == 0, result.stderr
    assert "from stdin" in result.stdout


def test_python_command_uses_plotnado_interpreter():
    import sys
    from pathlib import Path

    result = _run_cli("python", "-c", "import sys; print(sys.executable)")
    assert result.returncode == 0, result.stderr
    assert Path(result.stdout.strip()).resolve() == Path(sys.executable).resolve()


def test_python_command_passes_help_to_interpreter():
    result = _run_cli("python", "--help")
    assert result.returncode == 0
    assert "usage:" in result.stdout.lower()
    assert "Run Python with plotnado importable" not in result.stdout


@pytest.fixture
def jupyter_data_dir(tmp_path, monkeypatch):
    pytest.importorskip("ipykernel")
    monkeypatch.setenv("JUPYTER_DATA_DIR", str(tmp_path))
    return tmp_path


def _kernel_spec(data_dir, name):
    import json

    return json.loads((data_dir / "kernels" / name / "kernel.json").read_text())


def test_kernel_install_defaults(jupyter_data_dir):
    import sys

    from plotnado.cli.cli import get_version

    result = runner.invoke(app, ["kernel", "install"])
    assert result.exit_code == 0, result.output

    spec = _kernel_spec(jupyter_data_dir, "plotnado")
    assert spec["display_name"] == f"PlotNado {get_version()}"
    assert spec["argv"][0] == sys.executable


def test_kernel_install_custom_names(jupyter_data_dir):
    result = runner.invoke(
        app, ["kernel", "install", "--name", "pn-test", "--display-name", "My PlotNado"]
    )
    assert result.exit_code == 0, result.output
    assert _kernel_spec(jupyter_data_dir, "pn-test")["display_name"] == "My PlotNado"


def test_kernel_uninstall_removes_kernel(jupyter_data_dir):
    assert runner.invoke(app, ["kernel", "install", "--name", "pn-test"]).exit_code == 0
    result = runner.invoke(app, ["kernel", "uninstall", "--name", "pn-test"])
    assert result.exit_code == 0, result.output
    assert not (jupyter_data_dir / "kernels" / "pn-test").exists()


def test_kernel_uninstall_missing_kernel_fails(jupyter_data_dir):
    result = runner.invoke(app, ["kernel", "uninstall", "--name", "does-not-exist"])
    assert result.exit_code == 1
    assert "No Jupyter kernel named 'does-not-exist'" in result.output


@pytest.mark.parametrize("subcommand", ["install", "uninstall"])
def test_kernel_without_ipykernel_points_to_notebook_extra(subcommand, tmp_path, monkeypatch):
    import sys

    monkeypatch.setenv("JUPYTER_DATA_DIR", str(tmp_path))
    monkeypatch.setitem(sys.modules, "ipykernel", None)
    result = runner.invoke(app, ["kernel", subcommand])
    assert result.exit_code == 1
    assert "plotnado[notebook]" in result.output
    assert not (tmp_path / "kernels").exists()


def test_python_command_in_process_forwards_args_and_exit_code(monkeypatch):
    """Covers the command body in-process; the subprocess tests above check real behaviour."""
    import subprocess
    import sys

    calls = []

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 5)

    monkeypatch.setattr("plotnado.cli.env.subprocess.run", fake_run)
    result = runner.invoke(app, ["python", "script.py", "--flag", "-x"])
    assert result.exit_code == 5
    assert calls == [[sys.executable, "script.py", "--flag", "-x"]]


def test_get_version_falls_back_when_not_installed(monkeypatch):
    from importlib.metadata import PackageNotFoundError

    from plotnado.cli import cli as cli_module

    def missing(_name):
        raise PackageNotFoundError("plotnado")

    monkeypatch.setattr(cli_module, "_package_version", missing)
    assert cli_module.get_version() == "unknown"
    assert runner.invoke(app, ["--version"]).output.strip() == "plotnado unknown"


def test_main_entry_point_runs_app(monkeypatch, capsys):
    import sys

    from plotnado.cli.cli import get_version, main

    monkeypatch.setattr(sys, "argv", ["plotnado", "--version"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == f"plotnado {get_version()}"
