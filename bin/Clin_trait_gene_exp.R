#Usage: Rscript Clin_trait_gene_exp.R <FI:Exp.mat> <FI:in.clinical.data> <FI:gene.list> <DIR:output> 
args <- commandArgs(T)
library(ggplot2)
library(ggpubr)
library(reshape2)
expdata = read.table(args[1], header=TRUE, row.names = 1, sep = "\t")
traitData = read.table(args[2], header=TRUE, row.names = 1, sep = "\t")
gene_lst = read.table(args[3], header = F, sep = "\t")

Sample_IDs = colnames(expdata)
traitRows = match(Sample_IDs, rownames(traitData))
datTraits = traitData[traitRows,]
datTraits[datTraits==""] = NA
select_gene = match(gene_lst$V1, rownames(expdata))
select_gene = select_gene[!is.na(select_gene)]

MEs = as.data.frame(t(expdata[select_gene,]))
pdf(file = paste(args[4],"/Clin_trait_relat.pdf",sep=""), width = 8, height = 10)
for (i in 1:ncol(datTraits)){
  mat = as.data.frame(cbind(MEs, datTraits[,i]))
  colnames(mat)[ncol(mat)] = colnames(datTraits)[i]
  data = melt(mat,id.var=colnames(mat)[ncol(mat)])
  data = na.omit(data)
  data[,1] = as.character(data[,1])
  #   pdf(file = paste(args[3],"/Clinical/Module_clin_relationship.pdf",sep=""), width = 20, height = 12)
  if (length(levels(as.factor(as.character(data[,1])))) == 2){
    print(ggplot(data=data, aes(x=as.character(data[,1]),y=value))+
            geom_boxplot(aes(fill=data[,1]), outlier.size = 0.5, outlier.colour = "red", outlier.alpha = 0)+
            facet_wrap(~ variable, scales = "free")+
            labs(x= colnames(data)[1],y = "Gene expression")+
            guides(fill=guide_legend(title=colnames(data)[1]))+
            geom_jitter(position = 'jitter',size=0.5)+
            theme_bw() +
            stat_compare_means(label = "p.format", color = "red")+
            theme(strip.text=element_text(size=rel(1.5)),strip.background=element_rect(fill="#ABD3ED"),
                  axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1),
                  axis.title.x = element_text(size = 20),
                  axis.title.y = element_text(size = 20)))
  }else if(length(levels(as.factor(as.character(data[,1])))) > 2){
    print(ggplot(data=data, aes(x=as.character(data[,1]),y=value))+
            geom_boxplot(aes(fill=data[,1]), outlier.size = 0.5, outlier.colour = "red", outlier.alpha = 0)+
            facet_wrap(~ variable, scales="free")+
            labs(x= colnames(data)[1],y = "Gene expression")+
            guides(fill=guide_legend(title=colnames(data)[1]))+
            geom_jitter(position = 'jitter', size=0.2)+
            theme_bw() +
            stat_compare_means(label = "p.signif", method = "t.test", 
                               ref.group = levels(as.factor(as.character(data[,1])))[1], color = "red")+
            theme(strip.text=element_text(size=rel(1.5)),strip.background=element_rect(fill="#ABD3ED"),
                  axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1),
                  axis.title.x = element_text(size = 20),
                  axis.title.y = element_text(size = 20)))
  }
}
dev.off()



