#Usage: Rscript Vocano_plot.R deseq.xls out.prefix
args = commandArgs(T)
library("ggplot2")

res = read.table(args[1],header=T)
res <- na.omit(res)
tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj))

pdf(file = paste(args[2], ".pdf", sep=""))
par(mar = c(5, 6, 5, 5))
res$baseMean[res$baseMean>5000]=5000
res$baseMean[res$baseMean<10]=10
lfc = 1
pval = 0.05
nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))

plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col= alpha("black", 0))
points(tab[nosigGene, ], pch = 16, cex = res$baseMean/1000, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col= alpha("black", 0.1))
if (length(unique(signGenes_up)) > 1){
	points(tab[signGenes_up, ], pch = 16, cex = res$baseMean/1000, col = alpha("red", 0.4)) 
}
if (length(unique(signGenes_down)) > 1){
	points(tab[signGenes_down, ], pch = 16, cex = res$baseMean/1000, col = alpha("blue", 0.4)) 
}
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "orange", lty = 2) 
if (length(unique(signGenes_up)) > 1){
	text(tab[signGenes_up, ]$logFC, tab[signGenes_up, ]$negLogPval+3, row.names(res[signGenes_up,]), cex = 0.5, col = "red")
}
if (length(unique(signGenes_down)) > 1){
	text(tab[signGenes_down, ]$logFC, tab[signGenes_down, ]$negLogPval+3, row.names(res[signGenes_down,]), cex = 0.5, col = "blue")
}
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", 2, "fold"), paste("+", 2, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c(alpha("red", 0.6),alpha("blue", 0.6)))

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
	text(tab[signGenes_up, ]$logFC, tab[signGenes_up, ]$negLogPval+3, row.names(res[signGenes_up,]), cex = 0.5, col = "red")
}
if (length(unique(signGenes_down)) > 1){
	text(tab[signGenes_down, ]$logFC, tab[signGenes_down, ]$negLogPval+3, row.names(res[signGenes_down,]), cex = 0.5, col = "blue")
}
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", 2, "fold"), paste("+", 2, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))

plot(tab, pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col = alpha("black", 0))
points(tab[nosigGene, ], pch = 21, cex = res$baseMean/1000, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col = "black", bg = alpha("black", 0.1))
if (length(unique(signGenes_up)) > 1){ 
	points(tab[signGenes_up, ], pch = 21, cex = res$baseMean/1000, col = "black", bg = alpha("red", 0.3)) 
}
if (length(unique(signGenes_down)) > 1){
	points(tab[signGenes_down, ], pch = 21, cex = res$baseMean/1000, col = "black", bg = alpha("cornflowerblue", 0.3)) 
}
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "orange", lty = 2) 
if (length(unique(signGenes_up)) > 1){
	text(tab[signGenes_up, ]$logFC, tab[signGenes_up, ]$negLogPval+3, row.names(res[signGenes_up,]), cex = 0.5, col = "red")
}
if (length(unique(signGenes_down)) > 1){
	text(tab[signGenes_down, ]$logFC, tab[signGenes_down, ]$negLogPval+3, row.names(res[signGenes_down,]), cex = 0.5, col = "blue")
}
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", 2, "fold"), paste("+", 2, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))

#res$color = "black"
#res$color[signGenes_up] = "red"
#res$color[signGenes_down] = "blue"
#background<-theme_bw()
#mytheme<-theme(panel.border = element_blank(),axis.line = element_line(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#label<-labs(x="Fold Change",y="-log10 Pvalue")
#grids<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

#ggplot(res,aes(x=res$log2FoldChange,y=-log10(res$padj)))+
#  xlab("log2fold change")+ylab("-log10padj")+
#  geom_point(size=as.numeric(res$baseMean)/1000,color=as.factor(res$color),alpha=0.3)+
#  background+grids+mytheme+geom_hline(aes(yintercept=1),colour="red", linetype="dashed")+scale_fill_discrete()+
#  geom_vline(aes(xintercept=1),colour="red", linetype="dashed")+geom_vline(aes(xintercept=-1),colour="red", linetype="dashed")+
#  annotate("text",x = res$log2FoldChange[signGenes_up], y = -log10(res$padj[signGenes_up])+3,label = rownames(res)[signGenes_up],size=2.5,colour="red")+
#  annotate("text",x = res$log2FoldChange[signGenes_down], y = -log10(res$padj[signGenes_down])+3,label = rownames(res)[signGenes_down],size=2.5,colour="blue")+label

dev.off()
