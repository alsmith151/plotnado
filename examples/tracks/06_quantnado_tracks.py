from __future__ import annotations

from pathlib import Path

from plotnado import GenomicFigure


REPO_ROOT = Path(__file__).resolve().parents[2]
QUANTNADO_FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "quantnado"
ATAC_FIXTURE = QUANTNADO_FIXTURE_DIR / "atac_chr22.zarr"
RNA_FIXTURE = QUANTNADO_FIXTURE_DIR / "rna_chr22.zarr"
METH_FIXTURE = QUANTNADO_FIXTURE_DIR / "meth_chr22.zarr"
SNP_FIXTURE = QUANTNADO_FIXTURE_DIR / "snp_chr1.zarr"
SIGNAL_REGION = "chr22:10,520,211-10,520,711"
METHYLATION_REGION = "chr22:100-360"
VARIANT_REGION = "chr1:54,721-55,221"


def _require_fixture(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required example fixture is missing: {path}")
    return path


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    atac_path = _require_fixture(ATAC_FIXTURE)
    rna_path = _require_fixture(RNA_FIXTURE)
    meth_path = _require_fixture(METH_FIXTURE)
    snp_path = _require_fixture(SNP_FIXTURE)

    signal_fig = GenomicFigure(track_height=1.2)
    signal_fig.axis()
    signal_fig.quantnado_coverage(
        "atac",
        dataset_path=str(atac_path),
        title="ATAC coverage",
        color="#1f78b4",
        fill=True,
    )
    signal_fig.quantnado_stranded_coverage(
        "rna",
        dataset_path=str(rna_path),
        title="RNA stranded coverage",
        color="#d7301f",
        reverse_color="#4575b4",
    )

    signal_outfile = outdir / "track_quantnado_signal.png"
    signal_fig.save(signal_outfile, SIGNAL_REGION, dpi=220)
    print(f"Saved {signal_outfile}")

    methylation_fig = GenomicFigure(track_height=1.2)
    methylation_fig.axis()
    methylation_fig.quantnado_methylation(
        "meth",
        dataset_path=str(meth_path),
        title="Methylation percentage",
        color="#238b45",
        point_size=22,
    )

    methylation_outfile = outdir / "track_quantnado_methylation.png"
    methylation_fig.save(methylation_outfile, METHYLATION_REGION, dpi=220)
    print(f"Saved {methylation_outfile}")

    variant_fig = GenomicFigure(track_height=1.2)
    variant_fig.axis()
    variant_fig.quantnado_variant(
        "snp",
        dataset_path=str(snp_path),
        title="Variant allele balance",
        het_color="#756bb1",
        hom_alt_color="#d7301f",
        marker_size=38,
    )

    variant_outfile = outdir / "track_quantnado_variant.png"
    variant_fig.save(variant_outfile, VARIANT_REGION, dpi=220)
    print(f"Saved {variant_outfile}")


if __name__ == "__main__":
    main()