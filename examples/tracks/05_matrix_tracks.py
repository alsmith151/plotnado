from __future__ import annotations

from pathlib import Path

from plotnado import GenomicFigure


REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_DATA_DIR = REPO_ROOT / "tests" / "data"
COOLER_FIXTURE = TEST_DATA_DIR / "hg19.GM12878-MboI.matrix.2000kb.cool"
CAPCRUNCHER_FIXTURE = TEST_DATA_DIR / "capcruncher" / "SAMPLE-A_REP1.hdf5"
COOLER_REGION = "chr1:0-10,000,000"
CAPCRUNCHER_REGION = "chr14:69,878,303-69,879,094"


def _require_fixture(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required example fixture is missing: {path}")
    return path


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    cooler_path = _require_fixture(COOLER_FIXTURE)
    capcruncher_path = _require_fixture(CAPCRUNCHER_FIXTURE)

    cooler_fig = GenomicFigure(track_height=2.1)
    cooler_fig.axis()
    cooler_fig.cooler(
        str(cooler_path),
        title="GM12878 2 Mb cooler",
        balance=False,
        cmap="YlOrRd",
        max_value=60,
    )
    cooler_fig.cooler_average(
        [str(cooler_path), str(cooler_path)],
        title="Cooler average (same file twice)",
        balance=False,
        cmap="viridis",
        max_value=60,
    )

    cooler_outfile = outdir / "track_matrix_cooler_average.png"
    cooler_fig.save(cooler_outfile, COOLER_REGION, dpi=220)
    print(f"Saved {cooler_outfile}")

    capcruncher_fig = GenomicFigure(track_height=2.4)
    capcruncher_fig.axis()
    capcruncher_fig.capcruncher(
        str(capcruncher_path),
        title="CapCruncher viewpoint Slc25A37",
        viewpoint="Slc25A37",
        balance=False,
        cmap="magma",
    )

    capcruncher_outfile = outdir / "track_capcruncher.png"
    capcruncher_fig.save(capcruncher_outfile, CAPCRUNCHER_REGION, dpi=220)
    print(f"Saved {capcruncher_outfile}")


if __name__ == "__main__":
    main()