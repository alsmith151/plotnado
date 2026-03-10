"""Tests for generated Figure type stub."""

from __future__ import annotations

import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR_PATH = REPO_ROOT / "scripts" / "generate_figure_stub.py"
GENERATED_PATH = REPO_ROOT / "plotnado" / "figure.pyi"


def _load_generator_module():
    spec = importlib.util.spec_from_file_location("generate_figure_stub", GENERATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load generator module at {GENERATOR_PATH}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generated_figure_stub_is_up_to_date() -> None:
    module = _load_generator_module()
    expected_content = module.generate()
    actual_content = GENERATED_PATH.read_text()

    assert actual_content == expected_content, (
        "plotnado/figure.pyi is out of date. "
        "Run `python scripts/generate_figure_stub.py` from the repository root "
        "and commit the updated file."
    )
