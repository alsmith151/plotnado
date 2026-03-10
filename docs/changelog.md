# Changelog

## Unreleased

- Refactor to Pydantic-based track architecture.
- Added tiered BED/GTF/BigBed I/O utilities.
- Added `CoolerTrack`, `CapcruncherTrack`, `CoolerAverage`.
- Added `BigWigCollection` and `BigWigDiff`.
<<<<<<< HEAD
- Added QuantNado track integrations:
  - `QuantNadoCoverageTrack`
  - `QuantNadoStrandedCoverageTrack`
  - `QuantNadoMethylationTrack`
  - `QuantNadoVariantTrack`
  with new `GenomicFigure` helpers and aliases.
=======
>>>>>>> 2273572867bad7bb6edf1bf8f5ecff6cd4752d5b
- Added `GenomicFigure` conveniences (`plot_gene`, `plot_regions`, `extend`, TOML helpers).
- Breaking: renamed `Figure` to `GenomicFigure` in the public API.
- Added alias and migration documentation.
