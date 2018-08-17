message = "#=========================================================================================
#        USAGE: Rsctipt edgeR.R <Exp_count.mat> <Sample_info.txt> <outdir> [Design.txt]
#
#  DESCRIPTION: Analysis differential expression genes among different groups
#
#       INPUTS: <Exp_count.mat> : Reads count matrix (integral)
#               <Sample_info.txt> : Sample annotation file (Header: \"\tType\")
#               <outdir> : Output derectory. 
#               <Design.txt> : Compairision design (eg: A_vs_B\nA_vs_C) [Optional]
#
# REQUIREMENTS: [R Packages] : edgeR ggplot2
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 08/21/2017
#     REVISION: ---
#================================================================================================"

#if (length(args)<2){
#  stop(message)
#}
args = commandArgs(T)
library("edgeR")
library("ggplot2")

#args[1]="../1.Exp_Merge/Exp_count.round_off.mat"
#args[2]="../../tmp/Sample_info.txt"
countData <- read.table(args[1], header = TRUE, row.names = 1)
countData<-round(countData)
colData <- read.table(args[2], header = TRUE, row.names = 1)
fc = 2
lfc = 1
pval = 0.05
if (is.na(args[3])){
  args[3] = "."
}


if (is.na(args[4])){
  # List all possible pared comparision
  type_level <- levels(colData$Type)
  comb <- combn(type_level,2)
  
  #i = 1
  for (i in 1:length(comb[1,])){
    # Extract specific info for comparision
    colData1 <- subset(colData, Type == comb[1,i] | Type == comb[2,i])
    countData1 <- countData[,row.names(colData1)]
    # edgeR
    y <- DGEList(counts = countData1,group=colData1$Type)
    keep <- rowSums(cpm(y)>0) >= 3 #a CPM>1 in at least 3 samples
    y <- y[keep,]
    y <- calcNormFactors(y)
    design <- model.matrix(~as.factor(as.vector(colData1$Type)))
    rownames(design)<-colnames(y)
    y<-estimateDisp(y,design)
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,coef=2)
    FDR <- p.adjust(qlf$table$PValue, method="BH")
    summary(de<-decideTestsDGE(qlf))
    res <- data.frame(qlf$table, FDR)
    
    if (comb[1,i] > comb[2,i]){
      res$logFC <- -res$logFC
    }
    
    #  res$logFC <- -res$logFC  # -logFC
    res_up <- subset(res, logFC > lfc & FDR < pval)
    res_down <- subset(res, logFC < -lfc & FDR < pval)
    write.table(res, file = paste(args[3],"/", comb[2,i], "_vs_", comb[1,i], ".edgeR.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    write.table(res_up, file = paste(args[3],"/", comb[2,i], "_vs_", comb[1,i], ".edgeR.up_regulate.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    write.table(res_down, file = paste(args[3],"/", comb[2,i], "_vs_", comb[1,i], ".edgeR.down_regulate.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    res <- na.omit(res)
    # Vocano plot
    pdf(file = paste(args[3],"/", comb[2,i], "_vs_", comb[1,i], ".edgeR.VocanoPlot.pdf", sep=""))
    par(mar = c(5, 6, 5, 5))
    tab = data.frame(logFC = res$logFC, negLogPval = -log10(res$FDR))
    res$logCPM[res$logCPM>10]=10
    res$logCPM[res$logCPM<1]=1
    nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
    signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
    signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
    gap = max(tab$logFC)/50
    
    up_count = length(which(signGenes_up))
    down_count = length(which(signGenes_down))
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
    mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
    mtext(c(paste("-", fc, "fold"), paste("+", fc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
    mtext(c(comb[2,i], comb[1,i]), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line = 2)
    mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
    legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
    
    dev.off()
    
    # Reverse
    res$logFC <- -res$logFC  # -logFC
    res_up <- subset(res, logFC > lfc & FDR < pval)
    res_down <- subset(res, logFC < -lfc & FDR < pval)
    write.table(res, file = paste(args[3],"/", comb[1,i], "_vs_", comb[2,i], ".edgeR.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    write.table(res_up, file = paste(args[3],"/", comb[1,i], "_vs_", comb[2,i], ".edgeR.up_regulate.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    write.table(res_down, file = paste(args[3],"/", comb[1,i], "_vs_", comb[2,i], ".edgeR.down_regulate.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    res <- na.omit(res)
    # Vocano plot
    pdf(file = paste(args[3],"/", comb[1,i], "_vs_", comb[2,i], ".edgeR.VocanoPlot.pdf", sep=""))
    par(mar = c(5, 6, 5, 5))
    tab = data.frame(logFC = res$logFC, negLogPval = -log10(res$FDR))
    res$logCPM[res$logCPM>10]=10
    res$logCPM[res$logCPM<1]=1
    nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
    signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
    signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
    gap = max(tab$logFC)/50
    
    up_count = length(which(signGenes_up))
    down_count = length(which(signGenes_down))
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
    mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
    mtext(c(paste("-", fc, "fold"), paste("+", fc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
    mtext(c(comb[1,i], comb[2,i]), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line = 2)
    mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
    legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
    
    dev.off()
  }
}else{
  des = read.table(args[4], header = F)
  for(i in 1:nrow(des)){
    comb = strsplit(as.character(des[i, 1]), "_vs_")[[1]]
    colData1 <- subset(colData, Type == comb[1] | Type == comb[2])
    countData1 <- countData[,row.names(colData1)]
    # edgeR
    y <- DGEList(counts = countData1,group=colData1$Type)
    keep <- rowSums(cpm(y)>0) >= 3 #a CPM>1 in at least 3 samples
    y <- y[keep,]
    y <- calcNormFactors(y)
    design <- model.matrix(~as.factor(as.vector(colData1$Type)))
    rownames(design)<-colnames(y)
    y<-estimateDisp(y,design)
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,coef=2)
    FDR <- p.adjust(qlf$table$PValue, method="BH")
    summary(de<-decideTestsDGE(qlf))
    res <- data.frame(qlf$table, FDR)
    
    #Define compare derection
    if (comb[1] < comb[2]){
      res$logFC = -res$logFC
    }
    
    #  res$logFC <- -res$logFC  # -logFC
    res_up <- subset(res, logFC > lfc & FDR < pval)
    res_down <- subset(res, logFC < -lfc & FDR < pval)
    write.table(res, file = paste(args[3],"/", comb[1], "_vs_", comb[2], ".edgeR.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    write.table(res_up, file = paste(args[3],"/", comb[1], "_vs_", comb[2], ".edgeR.up_regulate.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    write.table(res_down, file = paste(args[3],"/", comb[1], "_vs_", comb[2], ".edgeR.down_regulate.xls",sep = ""),sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    res <- na.omit(res)
    # Vocano plot
    pdf(file = paste(args[3],"/", comb[1], "_vs_", comb[2], ".edgeR.VocanoPlot.pdf", sep=""))
    par(mar = c(5, 6, 5, 5))
    tab = data.frame(logFC = res$logFC, negLogPval = -log10(res$FDR))
    res$logCPM[res$logCPM>10]=10
    res$logCPM[res$logCPM<1]=1
    nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
    signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
    signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
    gap = max(tab$logFC)/50
    
    up_count = length(which(signGenes_up))
    down_count = length(which(signGenes_down))
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
    mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
    mtext(c(paste("-", fc, "fold"), paste("+", fc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
    mtext(c(comb[1], comb[2]), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line = 2)
    mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
    legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
    
    dev.off()
  }
}

