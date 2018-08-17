#Rscript CV_variant_gene_extract.R <Exp.logged.matrix> <output.predix>
args <- commandArgs(T)
require(DESeq2)
require(statmod)
M <- read.table(args[1],header=T)
M[1:5,1:5]
ed <- M
means <- rowMeans(ed)
ed <- ed[means > 0,] ### Remove unexpressed genes in all samples
vars <- apply(ed,1,var)
means <- rowMeans(ed)
cv2 <- vars/means^2
#minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], 0.95 ) )
minMeanForFit <- 1
useForFit <- means >= minMeanForFit
#useForFit <- means >= 10
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
#fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means ),cv2 )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
vfit <- a1/xg + a0
df <- ncol(ed) - 1
pdf(file=paste(args[2], "CV.pdf", sep="."))
plot(means,cv2,pch=16,col="#FFA50000",cex=0.8, ylim=c(0,0.2))
#smoothScatter(log(means)/log(10),log(cv2)/log(10),nrpoints = 0)
sum = 0
for (i in 1:length(means)){
  if (cv2[i] > (a1/means[i]+a0)*qchisq(0.975,df)/df){
    points(means[i],cv2[i],col="#FF000064",pch=19,cex=0.8)
    sum=sum+1
  }else{
    points(means[i],cv2[i],col="#80808064",pch=19,cex=0.8)
  }
}
lines(xg,vfit,col="black", lwd=3 )
lines(xg,vfit * qchisq(0.975,df)/df,lty=2,col="black")
lines(xg,vfit * qchisq(0.025,df)/df,lty=2,col="black")
vgen = M[which(cv2 > (a1/means+a0)*qchisq(0.975,df)/df),]
write.table(file=paste(args[2], ".variant_gene.txt", sep = ""),vgen,sep="\t",quote = FALSE,row.names = T,col.names = NA) 

dev.off()
