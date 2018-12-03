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
version="v0.0.1"
//=======================================================================================
// Nextflow Version check
if( !nextflow.version.matches('0.30+') ) {
    println print_yellow("This workflow requires Nextflow version 0.26 or greater -- You are running version ")+ print_red(nextflow.version)
    exit 1
}


params.help = null
if (params.help) {
    log.info ''
    log.info print_purple('------------------------------------------------------------------------')
    log.info "RNAseqPipe_SYSUCC:  v$version"
     log.info print_purple('------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ')
    log.info print_yellow('    The typical command for running the pipeline is as follows (we do not recommend users passing configuration parameters through command line, please modify the config.file instead):\n') +
            print_purple('       Nextflow run RNAseqPipe/main.nf \n') +

            print_yellow('    General arguments:             Input and output setting\n') +
            print_cyan('      --inputdir <path>            ') + print_green('Path to input data(optional), current path default\n') +
            print_cyan('      --reads <*_fq.gz>            ') + print_green('Filename pattern for pairing raw reads, e.g: *_{1,2}.fastq.gz for paired reads\n') +
            print_cyan('      --outdir <path>               ') + print_green('The output directory where the results will be saved(optional), current path is default\n') +
            '\n' +
            print_yellow('    Options:                         General options for run this pipeline\n') +
            print_cyan('      --design <file>               ') + print_green('A flat file stored the experimental design information ( required when perform differential expression analysis)\n') +
            print_cyan('      --compare <file>               ') + print_green('A flat file stored comparison information ( required when perform differential expression analysis, e.g )\n') +
            print_cyan('      --singleEnd                   ') + print_green('Reads type, True for single ended \n') +
            print_cyan('      --unstrand                    ') + print_green('RNA library construction strategy, specified for \'unstranded\' library \n') +
            '\n' +
            print_yellow('    References:                      If not specified in the configuration file or you wish to overwrite any of the references.\n') +
            print_cyan('      --fasta                       ') + print_green('Path to Fasta reference(required)\n') +
            print_cyan('      --gene_gtf                    ') + print_green('An annotation file from GENCODE database in GTF format (required)\n') +
           '\n' +
            print_yellow('    Other options:                   Specify the email and \n') +
             print_cyan('      --mail                         ') + print_green('email info for reporting status of your LncPipe execution  \n') +



            log.info '------------------------------------------------------------------------'
    log.info print_yellow('Contact information: zhaoqi@sysucc.org.cn')
    log.info print_yellow('Copyright (c) 2013-2018, Sun Yat-sen University Cancer Center.')
    log.info '------------------------------------------------------------------------'
    exit 0
}




// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
mail=params.email
// read file
datoolPath = file(params.dapath)
if( !datoolPath.exists() ) exit 1, print_red("DAtools  not found: ${params.dapath}")
gene_gtf = file(params.gtf)
if( !gene_gtf.exists() ) exit 1, print_red("GTF file not found: ${params.gtf}")
// star-rsem index
star_index =  params.star_index
gseapath = params.gseapath
gsea_pathway = params.gsea_pathway

//design file
if(params.designfile) {
    designfile = file(params.designfile)
    if( !designfile.exists() ) exit 1, print_red("Design file not found: ${params.designfile}")
}else{
    designfile=null
}
//compare.txt
if(params.comparefile){
    File comparefile = new File(params.comparefile)
    if( !comparefile.exists() ) exit 1, print_red("Compare file not found: ${params.comparefile}")
    compareLines = Channel.from(comparefile.readLines())
}else{
    compareLines=""
}
compareLines.into{compareLines_for_DE; compareLines_for_GSEA}



// Check parameters

//Checking parameters
log.info print_purple("You are running RNAseqPipe-SYSUCC with the following parameters:")
log.info print_purple("Checking parameters ...")
log.info print_yellow("=====================================")
log.info print_yellow("Fastq file extension:           ") + print_green(params.read)
log.info print_yellow("Single end :                    ") + print_green(params.singleEnd)
log.info print_yellow("Strand specific condition:      ") + print_green(params.strand)
log.info print_yellow("Output folder:                  ") + print_green(params.outdir)
log.info print_yellow("STAR index path:                ") + print_green(params.star_index)
log.info print_yellow("GTF path:                       ") + print_green(params.gtf)
log.info print_yellow("Design file  path:              ") + print_green(params.designfile)
log.info print_yellow("Compare file path:              ") + print_green(params.comparefile)
log.info print_yellow("=====================================")
log.info "\n"

/*
 Step : Fastqc by fastp
 */

reads = params.read

Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
        .ifEmpty {
    exit 1, print_red("Cannot find any reads matching: ${reads}\nPlz check your read string in nextflow.config file \n")
}
.set { reads_for_fastqc}


if(params.skip_qc){
    println print_yellow("Skip reads qc process, directly do mapping ")
    process RSEM_quantification  {

        tag { file_tag }
        label 'para'

        publishDir pattern: "*.bam",
                path: { params.outdir + "/Star_alignment" }, mode: 'move', overwrite: true

        input:
        set val(samplename), file(pair) from reads_for_fastqc

        output:
        file "${file_tag_new}.genes.results" into counting_file
        set val(file_tag_new), file ("${file_tag_new}.STAR.genome.bam") into bamfile_for_qualimap

        shell:
        println print_purple("Start analysis with RSEM  " + samplename)
        file_tag = samplename
        file_tag_new = file_tag


        println print_purple("Initial reads mapping of " + samplename + " performed by STAR in paired-end mode")
        """
             rsem-calculate-expression -p ${task.cpus} \
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
                path: { params.outdir + "/QC" }, mode: 'copy', overwrite: true

        input:
        set val(samplename), file(fastq_file) from reads_for_fastqc

        output:
        file "*.html" into fastqc_for_waiting,fastqc_for_multiqc
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

/*
 Step : Quantification  by star and RSEM
 */

    process RSEM_quantification_without_fastp {

        tag { file_tag }
        label 'para'

        publishDir pattern: "*.bam",
                path: { params.outdir + "/Star_alignment" }, mode: 'link', overwrite: true

        input:
        set val(samplename), file(pair) from readPairs_for_discovery
        file tempfiles from fastqc_for_waiting // just for waiting

        output:
        file "${file_tag_new}.genes.results" into counting_file
        set val(file_tag_new), file ("${file_tag_new}.STAR.genome.bam") into bamfile_for_qualimap
        shell:
        println print_purple("Start analysis with RSEM  " + samplename)
        file_tag = samplename
        file_tag_new = file_tag

    if(params.strand){
        println print_purple("Initial reads mapping of " + samplename + " performed by STAR in single-end mode")
        """
                 rsem-calculate-expression -p ${task.cpus} \
                    --no-bam-output --star \
                    -star-gzipped-read-file \
                    --star-output-genome-bam \
                    --estimate-rspd --time \
                    --strand-specific \
                    --paired-end  ${pair[0]} ${pair[1]} ${star_index} ${file_tag_new}
                    
        """
        }else{
        println print_purple("Initial reads mapping of " + samplename + " performed by STAR in paired-end mode")
        """
                     rsem-calculate-expression -p ${task.cpus} \
                        --no-bam-output --star \
                        -star-gzipped-read-file \
                        --star-output-genome-bam \
                        --estimate-rspd --time \
                        --paired-end  ${pair[0]} ${pair[1]} ${star_index} ${file_tag_new}
                        
        """
        }


    }
}

process run_qualimap{
    tag { file_tag }

    publishDir pattern: "*.pdf",
            path: { params.outdir + "/Qualimap" }, mode: 'copy', overwrite: true

    input:
    set val(samplename), file(bam_for_qualimap) from bamfile_for_qualimap
    file gene_gtf

    output:
    file "*" into qualimap_result_for_multiqc

    shell:
    file_tag = samplename
    file_tag_new = file_tag
    samplename = file_tag

    """
     qualimap rnaseq -bam ${bam_for_qualimap} -gtf ${gene_gtf} -s -outfile ${file_tag_new}_qualimap.pdf
    """

}


/*
  Merge Expression matrix
 */
process collapse_matrix{
    tag { file_tag }

    publishDir pattern: "*.matrix",
            path: { params.outdir + "/Express_matrix" }, mode: 'copy', overwrite: true

    input:
    file abundance_tsv_matrix from counting_file.collect()
    file gene_gtf

    output:
    file "${samplename}.count.matrix" into count_matrix
    file "forDE.count.matrix" into count_matrix_forDE
    file "${samplename}.tpm.matrix" into tpm_matrix_forDE
    file "${samplename}.fpkm.matrix" into fpkm_matrix_for_GSEA
    shell:
    file_tag = 'combine'
    file_tag_new = file_tag
    samplename = file_tag
    """
     java -jar ${datoolPath} -MM -mode RSEM ./ ${samplename}.count.matrix -count -gtf $gene_gtf
     java -jar ${datoolPath} -MM -mode RSEM ./ ${samplename}.tpm.matrix -tpm -gtf $gene_gtf
     java -jar ${datoolPath} -MM -mode RSEM ./ ${samplename}.fpkm.matrix -fpkm -gtf $gene_gtf
     perl  ${baseDir}/bin/Ensem2Symbol.pl ${samplename}.count.matrix > forDE.count.matrix


    """

}
/*
  Differential expression analysis
 */
if(params.designfile && params.comparefile){
    process Differential_Expression_analysis{
        tag {file_tag}

        publishDir pattern: "{*.mat,*.xls,*.pdf}",
                path: { params.outdir + "/DEG" }, mode: 'copy', overwrite: true

        input:

        file countMatrix from count_matrix_forDE
        file designfile
        val compare_str from compareLines_for_DE

        output:

        set val(comstr),file("${comstr}.deseq.xls") into DE_result
        file "*" into DE_result_out

        shell:
        comstr = compare_str
        file_tag = 'DESeq2: '+comstr
        file_tag_new = file_tag
        """
        ln -s ${baseDir}/bin/PCAplot.R .
        Rscript ${baseDir}/bin/DESeq2.R ${countMatrix} ${designfile} ./ ${comstr} 
        """

    }
/*
 Gene Set Enrichment Analysis
 */

    process GSEA_analysis{
        tag { file_tag }

        publishDir pattern: "*",
                path: { params.outdir + "/GSEA_analysis" }, mode: 'move', overwrite: true

        input:
        file fpkm_matrix from fpkm_matrix_for_GSEA
        file designfile
        val compare_str from compareLines_for_GSEA


        output:
        file "${compare_str}*" into gsea_out

        shell:
        file_tag = compare_str
        file_tag_new = file_tag
        samplename = file_tag

        """
        # generate GSEA rnk file 
         perl ${baseDir}/bin/get_preRankfile_for_GSEA.pl ${fpkm_matrix} ${designfile} ${compare_str} ${compare_str}.rnk
         java -cp ${gseapath} xtools.gsea.GseaPreranked \
          -gmx  ${gsea_pathway}\
          -norm meandiv -nperm 1000 \
          -rnk  ${compare_str}.rnk \
          -scoring_scheme weighted -rpt_label ${compare_str} \
          -create_svgs true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false \
          -out ./ -gui false
        """
    }

}

/*
MultiQC for data quality report
 */

process Run_MultiQC {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('*') from fastqc_for_multiqc.collect()
    file ('*') from qualimap_result_for_multiqc.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    when:
    !params.skip_multiqc

    script:
    """
	    multiqc --force --interactive .
	 """
}



/*
Working completed message
 */
workflow.onComplete {
    log.info print_green("=================================================")
    log.info print_green("Cheers! RNAseq Pipeline from SYSUCC run Complete!")
    log.info print_green("=================================================")
    //email information
    if (params.mail) {
        recipient = params.mail
        def subject = 'My RNAseq-SYSUCC execution'

        ['mail', '-s', subject, recipient].execute() <<
                """
            RNAseq-SYSUCC execution summary
            ---------------------------
            Your command line: ${workflow.commandLine}
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
        
            """
    }


}
workflow.onError {

    println print_yellow("Oops... Pipeline execution stopped with the following message: ")+print_red(workflow.errorMessage)
}