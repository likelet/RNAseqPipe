#!/usr/bin/env nextflow

/*
 * RNAseqPipe was implemented by Dr. Qi Zhao from Sun Yat-sen University Cancer Center.
 *
 *
 *   LncPipe is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *      See the GNU General Public License for more details.
 *
 *
 */

/*
 * RNAseqPipe: RNAseq analysis Pipeline for standard comparison
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>: design and implement the pipeline.
 */

// requirement:
// - fastp
// - STAR
// - RSEM
// - DAtools
//pre-defined functions for render command
//=======================================================================================
ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";


def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }

//Help information
// Nextflow  version
version="v0.2.31"
//=======================================================================================
// Nextflow Version check
if( !nextflow.version.matches('0.25+') ) {
    println print_yellow("This workflow requires Nextflow version 0.26 or greater -- You are running version ")+ print_red(nextflow.version)
    exit 1
}


// read file
datoolPath = file(params.dapath)

gene_gtf = file(params.gtf)
// star-rsem index
star_index =  params.star_index
//design file
designfile = file(params.designfile)
//compare.txt
File comparefile = new File(params.comp_file)
compareLines = Channel.from(comparefile.readLines())

/*
 Step : Fastqc by fastp
 */

reads = params.input_folder + params.fastq_ext

Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
        .ifEmpty {
    exit 1, print_red("Cannot find any reads matching: ${reads}\nPlz check your fastq_ext string in nextflow.config file \n")
}
.set { reads_for_fastqc}


if(params.skip_qc){
    println print_yellow("Skip reads qc process, directly do mapping ")
    process RSEM_quantification  {

        tag { file_tag }
        label 'para'

        publishDir pattern: "*.bam",
                path: { params.out_folder + "/Result/Star_alignment" }, mode: 'mv', overwrite: true

        input:
        set val(samplename), file(pair) from reads_for_fastqc

        output:
        file "${file_tag_new}.genes.results" into counting_file

        shell:
        println print_purple("Start analysis with RSEM  " + samplename)
        file_tag = samplename
        file_tag_new = file_tag
        star_threads = 40


        println print_purple("Initial reads mapping of " + samplename + " performed by STAR in single-end mode")
        """
             rsem-calculate-expression -p ${star_threads} \
                --no-bam-output --star \
                -star-gzipped-read-file \
                --star-output-genome-bam \
                --estimate-rspd --time \
                --paired-end  ${pair[0]} ${pair[1]} ${star_index} ${file_tag_new}
                
        """
    }
}else{
    process Run_FastP {

        tag { fastq_tag }
        label 'qc'

        publishDir pattern: "*.html",
                path: { params.out_folder + "/Result/QC" }, mode: 'copy', overwrite: true

        input:
        set val(samplename), file(fastq_file) from reads_for_fastqc

        output:
        file "*.html" into fastqc_for_waiting
        set val(fastq_tag), file('*qc.fq.gz')  into readPairs_for_discovery
        shell:
        fastq_tag = samplename
        if (params.singleEnd) {
            """
        fastp -i ${fastq_file[0]} -o ${samplename}.qc.gz -h ${samplename}_fastp.html
        """


        } else {
            """

        fastp -i ${fastq_file[0]}  -I ${fastq_file[1]} -o ${samplename}_1.qc.fq.gz  -O ${samplename}_2.qc.fq.gz -h ${samplename}_fastp.html
"""

        }
    }

    fastqc_for_waiting = fastqc_for_waiting.first()
/*
 Step : Quantification  by star and RSEM
 */

    process RSEM_quantification  {

        tag { file_tag }
        label 'para'

        publishDir pattern: "*.bam",
                path: { params.out_folder + "/Result/Star_alignment" }, mode: 'mv', overwrite: true

        input:
        set val(samplename), file(pair) from readPairs_for_discovery
        file tempfiles from fastqc_for_waiting // just for waiting

        output:
        file "${file_tag_new}.genes.results" into counting_file

        shell:
        println print_purple("Start analysis with RSEM  " + samplename)
        file_tag = samplename
        file_tag_new = file_tag
        star_threads = 40


        println print_purple("Initial reads mapping of " + samplename + " performed by STAR in single-end mode")
        """
             rsem-calculate-expression -p ${star_threads} \
                --no-bam-output --star \
                -star-gzipped-read-file \
                --star-output-genome-bam \
                --estimate-rspd --time \
                --paired-end  ${pair[0]} ${pair[1]} ${star_index} ${file_tag_new}
                
        """

    }
}


process collapse_matrix{
    tag { file_tag }

    publishDir pattern: "*.count.matix",
            path: { params.out_folder + "/Result/express_matrix" }, mode: 'copy', overwrite: true

    input:
    file abundance_tsv_matrix from counting_file.collect()
    file gene_gtf

    output:
    file "${samplename}.count.matrix" into count_matrix
    file "forDE.count.matrix" into count_matrix_forDE
    shell:
    file_tag = 'combine'
    file_tag_new = file_tag
    samplename = file_tag
    """
     java -jar ${datoolPath} -MM -mode RSEM ./ ${samplename}.count.matix -count
     perl  ${baseDir}/bin/Ensem2Symbol.pl ${samplename}.count.matix > forDE.count.matrix


    """


}

if(designfile!=null){
    process Differential_Expression_analysis{
        tag {file_tag}

        publishDir pattern: "{*.mat,*.xls,*.pdf}",
                path: { params.out_folder + "/Result/DEG" }, mode: 'copy', overwrite: true

        input:

        file countMatrix from count_matrix_forDE
        file designfile
        file comparestr from compareLines

        output:

        set val(comparestr),file("${comparestr}.deseq.xls") into DE_result

        shell:
        file_tag = 'DESeq2: '+comparestr
        file_tag_new = file_tag
        """
        ln -s ${baseDir}/bin/PCAplot.R 
        Rscript ${baseDir}/bin/DESeq2.R ${countMatrix} ${designfile} ./ ${comparestr} 
        """

    }

    process GO_kegg_enrichment_analysis{
        tag {file_tag}

        publishDir pattern: "",
                path: { params.out_folder + "/Result/GO_KEGG" }, mode: 'copy', overwrite: true

        input:

        file val(comparestr),file("${comparestr}.deseq.xls") from DE_result
        file designfile

        shell:
        file_tag = comparestr
        file_tag_new = file_tag
        """
        Rscript ${baseDir}/bin/GO_KEGG_Reactome_enrich_DEG.R ${comparestr}.deseq.xls
        """

    }

    process GSEA{

    }

}
