#Useage: Rscript Gene_TPM_dist.R Exp_TPM.all_0.TPM_dis.txt
args = commandArgs(T)
library(ggplot2)
library(reshape2)
#args[1] <- "Exp_TPM.all_0.TPM_dis.txt"
data <- read.delim(args[1], header=T)
data2 <- cbind(data[,12], data[,11], data[,10], data[,9], data[,8])
#data2 <- data[,8:12]
rownames(data2) <- data$Sample_ID
data3 <- as.data.frame(t(data2))
data3$type <- c("0<TPM<=0.1", "0.1<TPM<1", "1<TPM<=10", "10<TPM<=100", "TPM>100")
data3$mark <- c("a", "b", "c", "d", "e")
data4 <- as.data.frame(melt(data3, varnames = "type"))
name <- gsub(args[1], pattern = ".txt", replacement = "")
pdf(file = paste(name, "pdf", sep = "."), height = 6, width = 10)
ggplot(data=data4, aes(x=variable, y=value, fill=mark)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),axis.text.y = element_text(hjust = 0.5)) +
  xlab("Sample") + ylab("# of genes") +
  theme(axis.title=element_text(size=14)) +
  theme(panel.spacing = unit(0, "lines")) +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B"),
                    name = NULL,
                    breaks = c("e", "d", "c", "b", "a"),
                    labels = c("0<TPM<=0.1", "0.1<TPM<1", "1<TPM<=10", "10<TPM<=100", "TPM>100"))
dev.off()
