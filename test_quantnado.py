from plotnado import GenomicFigure, QuantNadoCoverageTrack
import quantnado as qn

ds_path = "/ceph/project/milne_group/cchahrou/processed_data/seqnado_output/multiomics/dataset/"
ds = qn.open_dataset(ds_path)

fig = GenomicFigure()
fig.autocolor()
fig.scalebar()
fig.genes("hg38")
fig.quantnado_coverage("ATAC-SEM-1", quantnado=ds, title="ATAC", title_color="black", normalize="rpkm", scaling_factor=0.4)
fig.quantnado_stranded_coverage("RNA-SEM-1", quantnado=ds, title="RNA", title_color="black")
fig.quantnado_methylation("TAPS-SEM", quantnado=ds, title="TAPS", title_color="black")
fig.quantnado_variant('gDNA-SEM', quantnado=ds, title="gDNA Variants", title_color="black")
plot = fig.plot_gene("GNAQ")

plot.savefig("test_quantnado.png", dpi=300)

