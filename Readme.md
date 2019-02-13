
# SYSUCC-RNAseqPipe

### Quick start  for reproductive analysis 

    nextflow run RNAseqPipe/main.nf -profile c2 --read "*_{1,2}.fq.gz" --designfile "design.file" --comparefile "compare.txt"

### Documentation
The SYSUCC-RNAseqPipe pipeline comes with documentation about the pipeline, found in the `docs/` directory:  

please find the information in [wiki](https://github.com/likelet/RNAseqPipe/wiki)  

1. [Installation and configuration](docs/Installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)



### Dependencies 
* Softwares 
    * [fastp](https://github.com/OpenGene/fastp)
    * [STAR](https://github.com/alexdobin/STAR)
    * [RSEM](https://deweylab.github.io/RSEM/)
    * [Qualimap](http://qualimap.bioinfo.cipf.es/)
    * [DAtools](https://github.com/likelet/DAtools)
    * [MultiQC](https://github.com/ewels/MultiQC)
    * [GSEA](http://software.broadinstitute.org/gsea/index.jsp)  
    * Several R packages for downstream analysis.

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
    
    suffix of your raw reads file. For example, `*_{1,2}.fq.gz` for paired end reads file `sampleA_1.fq.gz` and `sampleA_2.fq.gz `  
    
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
    
* `--without_replicate`   

    set `ture` if your have no biological replicate.       
    *note: for no replicate mode, the compare file should be directly specified as `SampleName_vs_SampleName` have just been trimmed by `read` suffix string *  
    
# Credits 
* Main author:
  * Qi Zhao ([@qi_likelet](https://github.com/likelet/))
* Contributors:
  * Xiaolong Zhang
  * Rucheng Diao
  * Kaiyu 
