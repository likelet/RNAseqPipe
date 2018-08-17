# useage: Rscript Fusion_circos.R <star-fusion.fusion_predictions.abridged.edit.tsv> <output.pdf>
args = commandArgs(T)
library(chimeraviz)
fusions <- import_starfusion(args[1], "hg38")
pdf(args[2])
plot_circle(fusions)
dev.off()
#fusion <- getFusionById(fusions, 3)
