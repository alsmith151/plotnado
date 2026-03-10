# Changelog

## Unreleased

- Refactor to Pydantic-based track architecture.
- Added tiered BED/GTF/BigBed I/O utilities.
- Added `CoolerTrack`, `CapcruncherTrack`, `CoolerAverage`.
- Added `BigWigCollection` and `BigWigDiff`.
- Added QuantNado track integrations:
  - `QuantNadoCoverageTrack`
  - `QuantNadoStrandedCoverageTrack`
  - `QuantNadoMethylationTrack`
  - `QuantNadoVariantTrack`
  with new `GenomicFigure` helpers and aliases.
- Added `GenomicFigure` conveniences (`plot_gene`, `plot_regions`, `extend`, TOML helpers).
- Breaking: renamed `Figure` to `GenomicFigure` in the public API.
- Added alias and migration documentation.
