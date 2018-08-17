message = "#=========================================================================================
#        USAGE: Rscript DESeq2.R <Exp_count.mat> <Sample_info.txt> <oudir> [Design.txt]
#
#  DESCRIPTION: Analysis differential expression genes among different groups
#               PS: Normalization first, then paired-wised DEG analysis
#
#       INPUTS: <Exp_count.mat> : Reads count matrix (integral)
#               <Sample_info.txt> : Sample annotation file (Header: \"\tType\")
#               <outdir> : Output derectory.
#               <Design.txt> : Compairision design (eg: A_vs_B\nA_vs_C) [Optional]
#
# REQUIREMENTS: [R Packages] : DESeq2 ggplot2
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn, modified by Qi zhao
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/16/2017
#     REVISION: ---
#================================================================================================"

#if (length(args)<2){
#  stop(message)
#}

args = commandArgs(T)
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("pheatmap"))
suppressMessages(library("ggrepel"))
suppressMessages(library("FactoMineR"))
suppressMessages(library("ggpubr"))

source("PCAplot.R")

countData <- read.table(args[1], header = TRUE, row.names = 1,check.names=FALSE)
countData <- round(countData)
colData <- read.table(args[2], header = TRUE, row.names = 1)
outdir <- args[3]
# keep the order consistency
countData<-countData[,row.names(colData)]
fc = 2
lfc = log2(fc)
pval = 0.05
#filter gene that expressed less than 3 samples
filter.gene=ncol(countData)/2

#vocanoPlot
vocano_plot = function(Sample_1 = "A", Sample_2 = "B", lfc = 0, pval = 0.05){
  par(mar = c(5, 6, 5, 5))
  tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj)) 
  res$baseMean[res$baseMean>5000]=5000
  res$baseMean[res$baseMean<10]=10
  nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
  signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
  signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
  up_count = length(which(signGenes_up))
  down_count = length(which(signGenes_down))
  gap = max(sort(tab[signGenes_up, ]$negLogPval, decreasing = T)[1], sort(tab[signGenes_down, ]$negLogPval, decreasing = T)[1])/50
  plot(tab, pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col = alpha("black", 0))
  points(tab[nosigGene, ], pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col = "black", bg = "grey")
  if (length(unique(signGenes_up)) > 1){
    points(tab[signGenes_up, ], pch = 21, col = "black", bg = "red") 
  }
  if (length(unique(signGenes_down)) > 1){
    points(tab[signGenes_down, ], pch = 21, col = "black", bg = "cornflowerblue") 
  }
  abline(h = -log10(pval), col = "green3", lty = 2) 
  abline(v = c(-lfc, lfc), col = "orange", lty = 2) 
  if (length(unique(signGenes_up)) > 1){
    text(tab[signGenes_up, ]$logFC, tab[signGenes_up, ]$negLogPval+gap, row.names(res[signGenes_up,]), cex = 0.5, col = "red")
  }
  if (length(unique(signGenes_down)) > 1){
    text(tab[signGenes_down, ]$logFC, tab[signGenes_down, ]$negLogPval+gap, row.names(res[signGenes_down,]), cex = 0.5, col = "blue")
  }
  mtext(paste("padj =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
  mtext(c(paste("-", fc, "fold"), paste("+", fc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
  mtext(c(Sample_1, Sample_2), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=2)
  mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
  legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
}

#MAplot
ReplotMAplot<-function(diff_express,Title){
  ggmaplot(diff_express, main = Title,
  fdr = 0.05, fc = 2, size = 0.4,
  palette = c("#B31B21", "#1465AC", "darkgray"),
  genenames = as.vector(diff_express$name),
  legend = "top", top = 20,
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = ggplot2::theme_minimal()) %>% ggexport(filename = paste(Title, "_MAplot.pdf",sep=""))
}


keep <- rowSums(countData > 0) >= filter.gene #a Count>0 in at least 3 samples
countData <- countData[keep,]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Type)
dds <- DESeq(dds)
norm_counts = counts(dds, normalized=TRUE)
write.table(norm_counts, file = paste(outdir, "All_sample_normalized_count.deseq.mat", sep = "/"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)

compare.str = args[4]
    comb = strsplit(as.character(compare.str), "_vs_")[[1]]
    res <- results(dds, contrast=c("Type",comb[1],comb[2]))
    res <- as.data.frame(res)
    write.table(res, file = paste(outdir,"/", comb[1], "_vs_", comb[2], ".deseq.xls",sep = ""),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
    res_up <- subset(res, log2FoldChange > lfc & padj < pval)
    write.table(res_up, file = paste(outdir,"/", comb[1], "_vs_", comb[2], ".deseq.up_regulate.xls",sep = ""),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
    res_down <- subset(res, log2FoldChange < -lfc & padj < pval)
    write.table(res_down, file = paste(outdir,"/", comb[1], "_vs_", comb[2], ".deseq.down_regulate.xls",sep = ""),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
    # Plot
    ReplotMAplot(res,paste(comb[1], "_vs_", comb[2],sep=""))
    pdf(file = paste(outdir,"/", comb[1], "_vs_", comb[2], ".deseq.Plot.pdf", sep=""))


    vocano_plot(Sample_1 = comb[1], Sample_2 = comb[2], lfc = lfc, pval = pval)
    deg <- subset(res, abs(log2FoldChange) > lfc & padj < pval)
    if (nrow(deg) > 2){
      colData1 <- subset(colData, Type == comb[1] | Type == comb[2])
      countData1 <- norm_counts[row.names(deg), row.names(colData1)]
      select <- order(rowMeans(countData1), decreasing = TRUE)
      countData1 = log2(countData1+1)
      pheatmap(countData1[select,order(colData1$Type)], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=colData1, scale = "row")
      pheatmap(countData1[select,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=colData1, scale = "row")
      pcaData <- as.data.frame(prcomp(countData1[select,])$rotation)

      pca_plot_text <- getPCAplot(pcaData,colData, isText=TRUE)
      print(pca_plot_text)
    }
    dev.off()





#PCA plot
