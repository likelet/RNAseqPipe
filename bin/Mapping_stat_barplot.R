#!/usr/bin/Rscript
#================================================================
#   filename:double.axis.R
#   author:Hanshi Xu
#   created time:2017.08.12
#   description:for mapping visualization
#   example: Rscript double.axis.R ./mapping.stat.txt ./mappling.stat.jpg
#   email:xhs.m.g@163.com
#================================================================

myarg = commandArgs(T)
if (length(myarg)!=1) stop("Rscript double.axis.R mapping.stat\n")
#myarg[1] = "/disk/zxl/projects/czh_duodenal_rna/4.MapQC/Map_stat/Map_stat.txt"
myda = read.table(myarg[1],header = T,sep="\t")
# color define start ------------------------------------------------------

  total_count = "#E64B35"
  uniq_count = "#4DBBD5"
  uniq_rate = "#01A087"
  muti_read = "#3C5C88"
  muti_count ="#F39B7F"
  first_axis = "black"
  secon_axis = "black"

# color define end --------------------------------------------------------

name <- gsub(myarg[1], pattern = ".txt", replacement = "")
pdf(file = paste(name, ".pdf", sep = ""), height = 6, width = 10)
par(mar=c(5,6,4,6))
myda$Uniq_map_rate = as.numeric(sub(pattern = "%",replacement = "",fixed = T,x=myda$Uniq_map_rate))
myda$Muti_map_count = as.numeric(sub(pattern = "%",replacement = "",fixed = T,x=myda$Muti_map_count))
conut = myda[,c(2,3,5)]
row.names(conut) = myda[,1]
rate = myda[c(4,6)]
row.names(rate) = myda[,1]
conut=t(as.matrix(conut))
conut=as.matrix(conut)
conut=conut/1000000
mconut=max(conut)
dconut=mconut*1
barplot(conut,beside = T,col=c(total_count,uniq_count,muti_read),xlim = c(0,4*length(conut[1,])),ylim=c(0,dconut),xaxt="n",yaxt="n")
axis(2, col=first_axis,col.axis=first_axis,las=1)
axis(4, col=secon_axis,col.axis=secon_axis,las=1,
  at=dconut/10 * (0:10),
  labels = (0:10)*10,
)

samnum=length(myda$Sample_ID)
xloc=4*((1:samnum) -1)+2.5
axis(1,labels = myda[,1],at=xloc,las=2)
par(new=T,mar=c(5,6,4,6))
plot(x=xloc,y=rate[,1],xlim = c(0,4*length(conut[1,])),xaxt="n" ,yaxt="n",bty="n",ylab = "",xlab = "", type="l",ylim=c(0,110),lwd=2,col=uniq_rate)

par(new=T,mar=c(5,6,4,6))
plot(x=xloc,y=rate[,2],xlim = c(0,4*length(conut[1,])),xaxt="n" ,yaxt="n",bty="n",ylab = "",xlab = "", type="l",ylim=c(0,110),lwd=2,col=muti_count)

mtext("Mapping ratio (%)",side=4,col=secon_axis,line=2.5) 
mtext("Reads count (M)",side=2,col=first_axis,line=2.5)
legend(
  "top",
  legend = colnames(myda[,-1]),
  text.col = c(total_count,uniq_count,uniq_rate,muti_read,muti_count),
  col=c(total_count,uniq_count,uniq_rate,muti_read,muti_count),
  pch = c(15,15,NA,15,NA),
  lty = c(NA,NA,1,NA,1),
  cex = 0.5,
  x.intersp = 0.7
)
dev.off()
