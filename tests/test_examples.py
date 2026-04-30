from __future__ import annotations

import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_RUNNER = REPO_ROOT / "examples" / "run_examples.py"


def test_example_runner_completes_for_local_examples() -> None:
    result = subprocess.run(
        [sys.executable, str(EXAMPLE_RUNNER)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        "examples/run_examples.py failed\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )