#Usage: Rscript Survival_gene_exp.R <FI:Exp.mat> <FI:in.clinical.data> <FI:gene.list> <DIR:output> 
args <- commandArgs(T)
library("survival")
library("survminer")
expdata = read.csv2(args[1], header=TRUE, row.names = 1, sep = "\t", fill = TRUE)
traitData = read.table(args[2], header=TRUE, row.names = 1, sep = "\t", fill = TRUE)
gene_lst = read.table(args[3], header = F, sep = "\t")

Sample_IDs = colnames(expdata)
traitRows = match(Sample_IDs, rownames(traitData))
datTraits = traitData[traitRows,]
datTraits[datTraits==""] = NA
select_gene = match(gene_lst$V1, rownames(expdata))
select_gene = select_gene[!is.na(select_gene)]

MEs = t(expdata[select_gene,])
pdf(file = paste(args[4],"/survival.pdf",sep=""), width = 20, height = 12)
for (i in seq(1,ncol(datTraits),2)){
  j = i + 1
  for(n in 1:ncol(MEs)){
#    mat = data.frame(datTraits[,i:j], "Add" = NA)
    mat = data.frame(datTraits[,i:j])
    mat = na.omit(mat)
    mat$Add = NA
    colnames(mat)[3] = colnames(MEs)[n]
    MEs1 = MEs[rownames(mat),]
    mat[,3][MEs1[,n] <= summary(as.numeric(MEs1[,n]))[2]] = "Low expression"
    mat[,3][MEs1[,n] >= summary(as.numeric(MEs1[,n]))[5]] = "High expression"
#    mat = na.omit(mat)
    fit<-survfit(Surv(mat[,1],mat[,2])~mat[,3],data=mat)
    ggsurv <- ggsurvplot(
      fit,                     # survfit object with calculated statistics.
      data = lung,             # data used to fit survival curves.
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      conf.int = TRUE,         # show confidence intervals for 
      # point estimates of survival curves.
      palette = c("#E7B800", "#2E9FDF"),
      #  xlim = c(0,500),         # present narrower X axis, but not affect
      # survival estimates.
      xlab = "Time in months",   # customize X axis label.
      #      break.time.by = 100,     # break X axis in time intervals by 500.
      ggtheme = theme_light(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T,# colour risk table text annotations.
      risk.table.height = 0.25, # the height of the risk table
      risk.table.y.text = FALSE,# show bars instead of names in text annotations
      # in legend of risk table.
      ncensor.plot = TRUE,      # plot the number of censored subjects at time t
      ncensor.plot.height = 0.25,
      #  conf.int.style = "step",  # customize style of confidence intervals
      #      surv.median.line = "hv",  # add the median survival pointer.
      legend.labs = levels(as.factor(mat[,3])) # change legend labels.
    )
    # Labels for Survival Curves (plot)
    ggsurv$plot <- ggsurv$plot + labs(
      title    = paste("Survival curves --- ", colnames(mat)[1], 
                       " vs gene ", colnames(mat)[3], sep =""),                     
      subtitle = "Based on Kaplan-Meier estimates", 
      y = colnames(mat)[1]
    )
    # Labels for Risk Table 
    ggsurv$table <- ggsurv$table + labs(
      title    = "Note the risk set sizes",          
      subtitle = "and remember about censoring."
    )
    # Labels for ncensor plot 
    ggsurv$ncensor.plot <- ggsurv$ncensor.plot + labs( 
      title    = "Number of censorings", 
      subtitle = "over the time."
    )
    print(ggsurv)
  }
}
dev.off()
  