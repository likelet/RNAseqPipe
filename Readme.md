
### Quick start  for reproductive analysis 

    nextflow -c nextflow.config run main.nf -resume -with-trace  -with-timeline timeline.html


### Dependencies 
* Softwares 
    * [fastp]
    * [STAR] 
    * [RSEM]
    * [Qualimap]
    * [DAtools](https://github.com/likelet/DAtools)
    * [GSEA]
* R packages 
 
        #bioconductor package
        dep.lib <- c("DESeq2","ReactomePA","CLusterProfile","org.Hs.eg.db","pathview","topGO")
        c.lib<-c("ggplot2","ggpubr","ggrepel","pheatmap","FactoMineR")



### Input file  

* `design.txt`  
sampleInfor presents the experimental design of your data set, it is just like a design file of `DESeq2` and `EdgeR` input.  

        Sample	Type
        P1003NA	N
        P1003TA	T
        P1162NA	N
        P1162TA	T
        P1408NA	N
        P1408TA	T
        P1527NA	N
        
* `compare.txt`
specify which group to compare in your differential expression analysis 
        
        T_vs_N
       
`T` and `N` are the identical strings as the `Type` column in `design.txt`.


### command parameters 



* `--reads`  
    
    suffix of your raw reads file. 
    
* `--designfile`  
    
    design file  
    
* `--comparefile`  
    
    compare file 
    
* `--gene_gtf`  
    
    gtf file for building your STAR index 

* `--singleEnd`  
    
    `true` when using a single End reads input, default `false` 

* `--strand`  
    
    `true` when using strand specific library , default `false` 
     
* `--skip_qc`   

    set `ture` if you are going to skip qc step 