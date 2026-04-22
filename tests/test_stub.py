"""Regression test: every public method on GenomicFigure must have a stub entry in figure.pyi.

This prevents the silent divergence that occurred when new methods were added to
figure.py but the generation script was not updated.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

from plotnado.figure import GenomicFigure

STUB_PATH = Path(__file__).resolve().parents[1] / "plotnado" / "figure.pyi"

# Methods intentionally omitted from the stub (private implementation details).
# _repr_html_ is included despite the underscore prefix — it's a Jupyter protocol method.
PRIVATE_SKIP = {
    "_apply_autocolor",
    "_apply_autoscale_groups",
    "_apply_extend",
    "_apply_theme_palette",
    "_aesthetic_explicitly_set",
    "_create_track_from_alias",
    "_font_family",
    "_is_meta_track",
    "_is_theme_palette_eligible",
    "_label_field_explicitly_set",
    "_resolve_gene_track",
    "_resolve_theme",
    "_set_auto_color",
    "_should_autocolor_track",
    "_track_registry",
}


def _stub_method_names() -> set[str]:
    tree = ast.parse(STUB_PATH.read_text())
    class_node = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.ClassDef) and node.name == "GenomicFigure"
    )
    return {
        node.name
        for node in ast.walk(class_node)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _public_methods() -> set[str]:
    include_despite_underscore = {"_repr_html_"}
    result = set()
    for name in dir(GenomicFigure):
        if name in PRIVATE_SKIP:
            continue
        if name.startswith("_") and name not in include_despite_underscore:
            continue
        raw = inspect.getattr_static(GenomicFigure, name)
        func = raw.__func__ if isinstance(raw, (classmethod, staticmethod)) else raw
        if callable(func):
            result.add(name)
    return result


def test_stub_covers_all_public_methods():
    stub_names = _stub_method_names()
    public_names = _public_methods()
    missing = public_names - stub_names
    assert not missing, (
        f"Public GenomicFigure methods not in figure.pyi: {sorted(missing)}\n"
        "Run `python scripts/generate_figure_stub.py` to regenerate the stub."
    )
