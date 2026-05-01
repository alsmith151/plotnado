from __future__ import annotations

import argparse
import importlib.util
import os
import runpy
import sys
import traceback
from pathlib import Path

CORE_SCRIPTS = [
    "basic_figure.py",
    "advanced_features.py",
    "quickstart/01_first_plot.py",
    "quickstart/02_aliases_and_options.py",
    "tracks/01_bigwig_styles.py",
    "tracks/02_bed_and_narrowpeak.py",
    "tracks/03_links_annotations.py",
    "recipes/01_autoscale_overlay_highlight.py",
    "recipes/02_theme_labels_toml.py",
    "recipes/03_gene_label_strategies.py",
]

OPTIONAL_SCRIPTS = [
    ("cooler", "tracks/05_matrix_tracks.py"),
    ("quantnado", "tracks/06_quantnado_tracks.py"),
]

REMOTE_SCRIPTS = [
    "tracks/04_bigwig_collection_and_diff.py",
]


def _module_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def scripts_to_run(include_remote: bool = False) -> list[str]:
    scripts = list(CORE_SCRIPTS)

    for module_name, rel_path in OPTIONAL_SCRIPTS:
        if _module_available(module_name):
            scripts.append(rel_path)
        else:
            print(f"Skipping {rel_path} because optional dependency '{module_name}' is not installed.")

    if include_remote:
        scripts.extend(REMOTE_SCRIPTS)
    else:
        print("Skipping remote examples. Pass --include-remote to run Blueprint-backed examples.")

    return scripts


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run PlotNado example scripts.")
    parser.add_argument(
        "--include-remote",
        action="store_true",
        help="Also run examples that stream remote public datasets.",
    )
    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parent
    repo_root = root.parent
    repo_root_str = str(repo_root)
    if repo_root_str not in sys.path:
        sys.path.insert(0, repo_root_str)
    failures: list[str] = []

    for rel_path in scripts_to_run(include_remote=args.include_remote):
        script = root / rel_path
        print(f"\n>>> Running {script.relative_to(root)}")
        previous_cwd = Path.cwd()
        try:
            os.chdir(repo_root)
            runpy.run_path(str(script), run_name="__main__")
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
            if code != 0:
                failures.append(rel_path)
        except Exception:
            traceback.print_exc()
            failures.append(rel_path)
        finally:
            os.chdir(previous_cwd)

    if failures:
        print("\nSome examples failed:")
        for rel_path in failures:
            print(f" - {rel_path}")
        raise SystemExit(1)

    print("\nAll examples completed successfully.")


if __name__ == "__main__":
    main()
