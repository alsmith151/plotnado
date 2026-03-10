from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPTS = [
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


def main() -> None:
    root = Path(__file__).resolve().parent
    failures: list[str] = []

    for rel_path in SCRIPTS:
        script = root / rel_path
        print(f"\\n>>> Running {script.relative_to(root)}")
        result = subprocess.run([sys.executable, str(script)], cwd=str(root.parent))
        if result.returncode != 0:
            failures.append(rel_path)

    if failures:
        print("\\nSome examples failed:")
        for rel_path in failures:
            print(f" - {rel_path}")
        raise SystemExit(1)

    print("\\nAll examples completed successfully.")


if __name__ == "__main__":
    main()
