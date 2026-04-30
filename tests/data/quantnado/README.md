Compact QuantNado fixture stores used by PlotNado's file-backed tests and examples.

- `atac_chr22.zarr`: created from QuantNado's packaged `MV411-CAT_MLL-N-1_subsample.bam` test BAM, restricted to `chr22`, sample name `atac`.
- `rna_chr22.zarr`: created from the same BAM as a stranded RNA store on `chr22`, sample name `rna`.
- `snp_chr1.zarr`: created from QuantNado's packaged `test.vcf.gz`, restricted to `chr1`, sample name `snp`.
- `meth_chr22.zarr`: small synthetic methylation store with the same per-sample layout expected by QuantNado's methylation reader, sample name `meth`.

These fixtures keep PlotNado's QuantNado coverage network-free while still exercising real on-disk store layouts.