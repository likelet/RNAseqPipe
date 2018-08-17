# useage: Rscript GO_KEGG_Reactome_enrich.R <DEG_result_updown.xls> <DEG_result.xls> [outdir]
args = commandArgs(T)
suppressMessages(library(clusterProfiler))
suppressMessages(library(ReactomePA))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(topGO))
suppressMessages(library(GOSemSim))
suppressMessages(library(magrittr))
suppressMessages(library(DOSE))
options(bitmapType='cairo')
options(stringsAsFactors = FALSE)

if (is.na(args[2])){
  outdir = "."
}else{
  outdir = args[2]
}

#args[1] <- "CSC_vs_CC.deseq.up_regulate.xls"
#args[2] <- "CSC_vs_CC.deseq.xls"

Name = basename(args[1])
Name = sub(".xls", "", Name)
Name = sub(".txt", "", Name)

use.pathview=FALSE

data <- read.table(args[1], header = T)
#data2 <- read.table(args[2], header = T)
deg_lst = which((data$log2FoldChange > 1 | data$log2FoldChange < -1) & data$padj < 0.05)
data_deg <- data[deg_lst,]

trans = bitr(rownames(data_deg), fromType="SYMBOL", toType=c("ENTREZID", "GENENAME"), OrgDb="org.Hs.eg.db")

geneList_1 <- data_deg$log2FoldChange
names(geneList_1) <- rownames(data_deg)
geneList_2 <- data$log2FoldChange
names(geneList_2) <- rownames(data)


# GO analysis (BP)

GO_analysis = function(genes = rownames(data_deg), type = "BP", pval = 0.01, qval = 0.05, outdir = ".", Name = NULL){
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = type,
                  pvalueCutoff  = pval,
                  qvalueCutoff  = qval)
  if (nrow(ego) > 0){
    write.table(summary(ego), file = paste(outdir, "/", Name, ".GO_", type, "_enrich.xls",sep = ""),
                sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    str_length = max(nchar(summary(ego)$Description))
    str_height = nrow(ego)
    if (str_height >20){
      pdf(paste(outdir,"/", Name, ".GO_", type, "_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 4)
    }else if(str_height < 6){
      pdf(paste(outdir,"/", Name, ".GO_", type, "_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 2.5)
    }else{
      pdf(paste(outdir,"/", Name, ".GO_", type, "_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 0.4*str_height)
    }
    #    print(barplot(ego, showCategory = nrow(ego)))
    print(barplot(ego, showCategory = 20))
    #    print(dotplot(ego, showCategory = nrow(ego)))
    print(dotplot(ego, showCategory = 20))
    dev.off()
    pdf(paste(outdir,"/", Name, ".GO_", type, "_enrich_plot2.pdf",sep = ""), width = 25, height = 25)
    emapplot(ego)
    cnetplot(ego, categorySize="pvalue", foldChange=data_deg$foldChange)
    plotGOgraph(ego)
    dev.off()
  }
}

print("============= GO BP ================")
GO_analysis(genes = rownames(data_deg), type = "BP", outdir = outdir, Name = Name)
print("============= GO MF ================")
GO_analysis(genes = rownames(data_deg), type = "MF", outdir = outdir, Name = Name)
print("============= GO CC ================")
GO_analysis(genes = rownames(data_deg), type = "CC", outdir = outdir, Name = Name)


#####################

#ego <- enrichGO(gene          = rownames(data),
#                OrgDb         = org.Hs.eg.db,
#                keytype       = 'SYMBOL',
#                ont           = "all",
#                pvalueCutoff  = 0.01,
#                qvalueCutoff  = 0.05)
#print(barplot(ego, split = "ONTOLOGY", showCategory = 20) + facet_grid(ONTOLOGY~.,scales = "free"))
####################

# Reactome analysis
print("============= Reactome ================")
re <- enrichPathway(gene         = trans$ENTREZID,
                    organism     = 'human',
                    pvalueCutoff = 0.05)
re <- setReadable(re, OrgDb = org.Hs.eg.db, keytype = "ENTREZID")
if (nrow(re) > 0){
  write.table(summary(re), file = paste(outdir,"/", Name, ".Reactome_enrich.xls",sep = ""),
              sep = "\t", quote = FALSE, row.names = T, col.names = NA)
  str_length = max(nchar(summary(re)$Description))
  str_height = nrow(re)
  if (str_height > 20){
    pdf(paste(outdir,"/", Name, ".Reactome_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 4)
  }else if(str_height < 6){
    pdf(paste(outdir,"/", Name, ".Reactome_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 2.5)
  }else{
    pdf(paste(outdir,"/", Name, ".Reactome_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 0.4*str_height)
  }
  #  print(barplot(re, showCategory = nrow(re)))
  #  print(dotplot(re, showCategory = nrow(re)))
  print(barplot(re, showCategory = 20))
  print(dotplot(re, showCategory = 20))
  pdf(paste(outdir,"/", Name, ".Reactome_enrich_plot2.pdf",sep = ""), width = 25, height = 25)
  emapplot(re)
  dev.off()
}

# KEGG analysis
print("============= KEGG ================")
kk <- enrichKEGG(gene         = trans$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keytype = "ENTREZID")
if (nrow(kk) > 0){
  write.table(summary(kk), file = paste(outdir,"/", Name, ".KEGG_enrich.xls",sep = ""),
              sep = "\t", quote = FALSE, row.names = T, col.names = NA)
  str_length = max(nchar(summary(kk)$Description))
  str_height = nrow(kk)
  if (str_height > 20){
    pdf(paste(outdir,"/", Name, ".KEGG_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 4)
  }else if(str_height < 6){
    pdf(paste(outdir,"/", Name, ".KEGG_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 2.5)
  }else{
    pdf(paste(outdir,"/", Name, ".KEGG_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 0.4*str_height)
  }
  #  print(barplot(kk, showCategory = nrow(kk)))
  #  print(dotplot(kk, showCategory = nrow(kk)))
  print(barplot(kk, showCategory = 20))
  print(dotplot(kk, showCategory = 20))
  dev.off()
  pdf(paste(outdir,"/", Name, ".KEGG_enrich_plot2.pdf",sep = ""), width = 25, height = 25)
  emapplot(kk)
  dev.off()

    #use pathview
    if(use.pathview){
        require(pathview)
        for (i in 1:nrow(kk)){
    #browseKEGG(kk, 'hsa04392')
            pathway.id = summary(kk)$ID[i]
            pathview(gene.data  = geneList_1,
            pathway.id = pathway.id,
            kegg.dir = outdir,
            gene.idtype ="SYMBOL",
            out.suffix = paste(Name, "_enrich", sep = ""),
            species    = "hsa",
            kegg.native = TRUE,
            limit      = list(gene=max(abs(geneList_1)), cpd=1))

            pathview(gene.data  = geneList_2,
            pathway.id = pathway.id,
            kegg.dir = outdir,
            gene.idtype ="SYMBOL",
            out.suffix = paste(Name, "_enrich_all", sep = ""),
            species    = "hsa",
            kegg.native = TRUE,
            limit      = list(gene=max(abs(geneList_2)), cpd=1))

            pathview(gene.data   = geneList_1,
            pathway.id  = pathway.id,
            kegg.dir    = outdir,
            gene.idtype ="SYMBOL",
            out.suffix  = paste(Name, "_enrich", sep = ""),
            species     = "hsa",
            kegg.native = F,
            same.layer  = F,
            limit       = list(gene=max(abs(geneList_1)), cpd=1))

            pathview(gene.data   = geneList_2,
            pathway.id  = pathway.id,
            kegg.dir    = outdir,
            gene.idtype ="SYMBOL",
            out.suffix  = paste(Name, "_enrich_all", sep = ""),
            species     = "hsa",
            kegg.native = F,
            same.layer  = F,
            limit       = list(gene=max(abs(geneList_2)), cpd=1))

            system(paste("mv ", pathway.id, "* ", outdir, sep = ""))
    }
    }

}



