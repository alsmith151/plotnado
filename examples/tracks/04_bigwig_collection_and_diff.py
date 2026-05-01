from __future__ import annotations

import tempfile
from pathlib import Path

from plotnado import GenomicFigure
from plotnado.tracks import GenomicRegion
from plotnado.tracks.bigwig import BigWigTrack

import pybigtools


BLUEPRINT_BIGWIGS = {
    "T12-15 plasma cell RNA": "http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/tonsil/T12-15/plasma_cell/RNA-Seq/IDIBAPS/T12-15-PC.signal.star_grape2_crg.GRCh38.20150815.bw",
    "T12-16 plasma cell RNA": "http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/tonsil/T12-16/plasma_cell/RNA-Seq/IDIBAPS/T12-16-PC.signal.star_grape2_crg.GRCh38.20150815.bw",
}
REGION_CHROM = "chr1"
REGION_START = 155_190_000
REGION_END = 155_220_000
REGION = "chr1:155,190,000-155,220,000"


def _cache_bigwig_crop(label: str, url: str, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{label.lower().replace(' ', '_').replace('-', '_')}.bw"
    if cache_path.exists():
        return cache_path

    region = GenomicRegion(chromosome=REGION_CHROM, start=REGION_START, end=REGION_END)
    data = BigWigTrack(data=url).fetch_data(region)
    if data.empty:
        raise ValueError(f"No Blueprint signal fetched for {label} in {REGION}")

    records = [
        (REGION_CHROM, int(row.start), int(row.end), float(row.value))
        for row in data.itertuples(index=False)
    ]
    pybigtools.open(str(cache_path), "w").write({REGION_CHROM: REGION_END}, records)
    return cache_path


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    labels = list(BLUEPRINT_BIGWIGS.keys())
    with tempfile.TemporaryDirectory(prefix="plotnado-blueprint-") as tmpdir:
        cache_dir = Path(tmpdir)
        files = [
            str(_cache_bigwig_crop(label, url, cache_dir))
            for label, url in BLUEPRINT_BIGWIGS.items()
        ]

        fig = GenomicFigure(track_height=1.15)
        fig.axis()
        fig.scalebar()
        fig.bigwig_collection(
            files,
            title="Blueprint plasma-cell RNA collection",
            labels=labels,
            colors=["#d55e00", "#0072b2"],
            alpha=0.75,
        )
        fig.bigwig_diff(
            files[0],
            files[1],
            title="T12-15 minus T12-16",
            positive_color="#b22222",
            negative_color="#1f5aa6",
            bar_alpha=0.5,
        )

        outfile = outdir / "track_bigwig_collection_diff.png"
        fig.save(outfile, REGION, dpi=220)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    main()