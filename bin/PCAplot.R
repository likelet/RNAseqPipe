## PCA analysis plot

suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("FactoMineR"))

#define theme

#convert userinput data and condition list for PCA analysis


#get ggplot2 output result
getPCAplot <- function(data,conditionlist,isText=FALSE){

    plotDefaultTheme <- theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, color = "black",size = 15),
    axis.text.y = element_text(angle = 00, color = "black",size = 15),
    axis.title = element_text(face = "bold", color = "black", size = 20),
    legend.title = element_text(face = "bold", color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20)
    )
    dataForPCAinitialize<-function(data,colData){
        data<-t(data)
        data<-data.frame(condition=colData,data)
        return(data)
    }
    if(isText){
        a<-dataForPCAinitialize(data,conditionlist)
        pca <-PCA(a[,2:ncol(a)], scale.unit=T, graph=F)
        PC1 <- pca$ind$coord[,1]
        PC2 <- pca$ind$coord[,2]
        plotdata <- data.frame(Condition=a[,1],PC1,PC2)
        plotdata$Condition <- factor(plotdata$Condition)
        plot <- ggplot(plotdata, aes(PC1, PC2),environment = environment()) +
            geom_point(aes(colour = Condition,shape = Condition),size = 5) +
            geom_text_repel(aes(label=rownames(plotdata)), size=5, hjust=0.5, vjust=-0.5)+
            plotDefaultTheme+
            scale_fill_brewer(palette="Spectral")
        return(plot)
    }else{
        a<-dataForPCAinitialize(data,conditionlist)
        pca <-PCA(a[,2:ncol(a)], scale.unit=T, graph=F)
        PC1 <- pca$ind$coord[,1]
        PC2 <- pca$ind$coord[,2]
        plotdata <- data.frame(Condition=a[,1],PC1,PC2)
        plotdata$Condition <- factor(plotdata$Condition)
        plot <- ggplot(plotdata, aes(PC1, PC2),environment = environment()) +
            geom_point(aes(colour = Condition,shape = Condition),size = 5) +
            plotDefaultTheme+
            scale_fill_brewer(palette="Spectral")
        return(plot)

    }


}