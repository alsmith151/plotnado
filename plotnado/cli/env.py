"""Commands that expose plotnado's own Python environment.

These make the Python API usable after an isolated install such as
``uv tool install plotnado``, where ``import plotnado`` is not available
from the user's system Python.
"""

import subprocess
import sys

import typer

from .cli import app, get_version

KERNEL_NAME = "plotnado"

kernel_app = typer.Typer(
    help="Register plotnado's environment as a Jupyter kernel.",
    no_args_is_help=True,
)
app.add_typer(kernel_app, name="kernel")


@app.command(
    "python",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
    add_help_option=False,
)
def python_command(ctx: typer.Context) -> None:
    """Run Python with plotnado importable (scripts, -c, -m, or a REPL).

    All arguments are passed to the Python interpreter, e.g.

        plotnado python my_plot.py
        plotnado python -c "import plotnado; print(plotnado.__version__)"
    """
    result = subprocess.run([sys.executable, *ctx.args])
    raise typer.Exit(result.returncode)


def _require_ipykernel():
    try:
        from ipykernel import kernelspec
    except ImportError:
        typer.echo(
            "ipykernel is not installed in plotnado's environment.\n"
            "Reinstall with the notebook extra, e.g.\n"
            '  uv tool install --force "plotnado[notebook]"\n'
            '  pip install "plotnado[notebook]"',
            err=True,
        )
        raise typer.Exit(1)
    return kernelspec


@kernel_app.command("install")
def kernel_install(
    name: str = typer.Option(KERNEL_NAME, help="Kernel name used by Jupyter."),
    display_name: str = typer.Option(
        None, help="Name shown in Jupyter's kernel picker [default: PlotNado <version>]."
    ),
) -> None:
    """Install a Jupyter kernel that runs in plotnado's environment.

    Afterwards, pick the kernel in Jupyter Lab, Notebook or VS Code to
    `import plotnado` in notebooks.
    """
    kernelspec = _require_ipykernel()
    display = display_name or f"PlotNado {get_version()}"
    path = kernelspec.install(user=True, kernel_name=name, display_name=display)
    typer.echo(f"Installed Jupyter kernel '{display}' ({name}) at {path}")


@kernel_app.command("uninstall")
def kernel_uninstall(
    name: str = typer.Option(KERNEL_NAME, help="Kernel name to remove."),
) -> None:
    """Remove the plotnado Jupyter kernel."""
    _require_ipykernel()
    from jupyter_client.kernelspec import KernelSpecManager, NoSuchKernel

    try:
        KernelSpecManager().remove_kernel_spec(name)
    except NoSuchKernel:
        typer.echo(f"No Jupyter kernel named '{name}' is installed.", err=True)
        raise typer.Exit(1)
    typer.echo(f"Removed Jupyter kernel '{name}'")
