
### Dependencies 
* Softwares 
    * [fastp]
    * [STAR] 
    * [RSEM]
    * [DAtools](https://github.com/likelet/DAtools)
    * [gseapy](https://pypi.org/project/gseapy/) 
* R packages 
 
        #bioconductor package
        dep.lib <- c("DESeq2","ReactomePA","CLusterProfile","org.Hs.eg.db","pathview","topGO")
        c.lib<-c("ggplot2","ggpubr","ggrepel","pheatmap","FactoMineR")



### Input file  

* `sampleInfor.txt`  
sampleInfor presents the experimental design of your data set, it is just like a design file of `DESeq2` and `EdgeR` input.  

        	Type
        P1003NA	N
        P1003TA	T
        P1162NA	N
        P1162TA	T
        P1408NA	N
        P1408TA	T
        P1527NA	N
        
* `compare.txt`
specify which group to compare in your differential expression analysis 
        
        A_vs_B
        B_vd_A
`A` and `B` are the identical strings as the `Type` column in `sampleInfor.txt`.


### parameters 

* `skip_qc`   

    set `ture` if you are going to skip qc step 

* `fastq_ext`  
    
    suffix of your raw reads file. 

