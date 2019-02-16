#!/usr/bin/env nextflow

/*
 * RNAseqPipe was implemented by Dr. Qi Zhao from Sun Yat-sen University Cancer Center.
 *
 *
 *   RNAseqPipe is free software: you can redistribute it and/or modify
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

//Help information
// Nextflow  version
version="0.0.1"
//=======================================================================================
// Nextflow Version check
if( !nextflow.version.matches('0.30+') ) {
    println LikelikeUtils.print_yellow("This workflow requires Nextflow version 0.26 or greater -- You are running version ")+ LikelikeUtils.print_red(nextflow.version)
    exit 1
}


params.help = null
if (params.help) {
    this.helpMessage()
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
gene_gtf = file(params.gene_gtf)
if( !gene_gtf.exists() ) exit 1, print_red("GTF file not found: ${params.gene_gtf}")
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
compareLines.into{compareLines_for_DE; compareLines_for_GSEA;compareLines_for_DE_without_REP}


//Checking parameters
minimalInformationMessage()


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
    process RSEM_quantification_without_fastp  {

        tag { file_tag }
        label 'para'

        publishDir pattern: "*.bam",
                path: { params.outdir + "/Star_alignment" }, mode: 'link', overwrite: true

        input:
        set val(samplename), file(pair) from reads_for_fastqc

        output:
        file "${file_tag_new}.genes.results" into counting_file,couting_file_DE_without_Rep
        set val(file_tag_new), file ("${file_tag_new}.STAR.genome.bam") into bamfile_for_qualimap

        shell:
        println print_purple("Start analysis with RSEM  " + samplename)
        file_tag = samplename
        file_tag_new = file_tag


        if(params.singleEnd){
            println print_purple("Initial reads mapping of " + samplename + " performed by STAR in single-end  mode")
             if(params.strand){
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in Strand specific  mode")
                """
                        rsem-calculate-expression -p ${task.cpus} \
                            --no-bam-output --star \
                            -star-gzipped-read-file \
                            --star-output-genome-bam \
                            --estimate-rspd --time \
                            --strand-specific \
                            ${pair[0]} ${star_index} ${file_tag_new}
                            
                """
            }else{
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in non strand  mode")
                """
                    rsem-calculate-expression -p ${task.cpus} \
                        --no-bam-output --star \
                        -star-gzipped-read-file \
                        --star-output-genome-bam \
                        --estimate-rspd --time \
                        ${pair[0]} ${star_index} ${file_tag_new}
                        
            """
            }
        } else{
            println print_purple("Initial reads mapping of " + samplename + " performed by STAR in paired-end  mode")
                
            if(params.strand){
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in Strand specific  mode")
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
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in non strand  mode")
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
}else{


/*
 Step : Quality control
 */
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
    process RSEM_quantification {

        tag { file_tag }
        label 'para'

        publishDir pattern: "*.bam",
                path: { params.outdir + "/Star_alignment" }, mode: 'link', overwrite: true

        input:
        set val(samplename), file(pair) from readPairs_for_discovery
        file tempfiles from fastqc_for_waiting // just for waiting

        output:
        file "${file_tag_new}.genes.results" into counting_file,couting_file_DE_without_Rep
        set val(file_tag_new), file ("${file_tag_new}.STAR.genome.bam") into bamfile_for_qualimap
        shell:
        println print_purple("Start analysis with RSEM  " + samplename)
        file_tag = samplename
        file_tag_new = file_tag

    if(params.singleEnd){
            println print_purple("Initial reads mapping of " + samplename + " performed by STAR in single-end  mode")
             if(params.strand){
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in Strand specific  mode")
                """
                        rsem-calculate-expression -p ${task.cpus} \
                            --no-bam-output --star \
                            -star-gzipped-read-file \
                            --star-output-genome-bam \
                            --estimate-rspd --time \
                            --strand-specific \
                            ${pair[0]} ${star_index} ${file_tag_new}
                            
                """
            }else{
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in non strand  mode")
                """
                    rsem-calculate-expression -p ${task.cpus} \
                        --no-bam-output --star \
                        -star-gzipped-read-file \
                        --star-output-genome-bam \
                        --estimate-rspd --time \
                        ${pair[0]} ${star_index} ${file_tag_new}
                        
            """
            }
        } else{
            println print_purple("Initial reads mapping of " + samplename + " performed by STAR in paired-end  mode")
                
            if(params.strand){
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in Strand specific  mode")
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
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in non strand  mode")
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


    when:
    !params.skip_qualimap


    shell:
    file_tag = samplename
    file_tag_new = file_tag
    samplename = file_tag
    parsed_mem = task.memory.getGiga()+"G"
    """
     qualimap rnaseq -bam ${bam_for_qualimap} -gtf ${gene_gtf} -s -outfile ${file_tag_new}_qualimap.pdf --java-mem-size=${parsed_mem}
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
    file "${samplename}.count.matrix" into count_matrix, count_for_IDEA
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
 generate IDEA input file 
*/

process IDEA_Input_generate{
publishDir pattern: "*.csv",
            path: { params.outdir + "/IDEA_Input" }, mode: 'move', overwrite: true

    input:
    file count_matrix from count_for_IDEA
    file designfile

    output:
    file "*"

    shell:
    
    '''
     sed 's/\t/,/g' !{count_matrix} > count.matrix.csv
     sed 's/\t/,/g' !{designfile} > design.csv
    '''

}
/*
  Differential expression analysis && GSEA
 */
if( params.designfile && params.comparefile ){
// DE
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

        when:
        !params.skip_gsea

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
 Differential expression analysis without replicates
 */
if(params.comparefile && params.without_replicate) {
    process Differential_Expression_analysis_without_RF {
        tag { file_tag }

        publishDir pattern: "{*}",
                path: { params.outdir + "/DEG_without_Rep" }, mode: 'copy', overwrite: true

        input:

        file abundance_tsv_matrix from couting_file_DE_without_Rep.collect()
        file gene_gtf
        val compare_str from compareLines_for_DE_without_REP

        output:

        file "*"

        shell:
        comstr = compare_str
        file_tag = 'PoissonTest: ' + comstr
        file_tag_new = file_tag
        comstr_a = comstr.split("_vs_")
        """
        java -jar ${datoolPath} -RNAseq -mode poissonDE \
                                ${comstr_a[0]}.genes.results  ${comstr_a[1]}.genes.results ${comstr} \
                                -gtf ${gene_gtf}
        """

    }
}


/*
MultiQC for data quality report
 */

process Run_MultiQC {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:

    file ('*') from qualimap_result_for_multiqc.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    when:
    !params.skip_multiqc & !params.skip_qualimap

    script:
    """
	    multiqc --force --interactive .
	 """
}



/*
Working completed message
 */
workflow.onComplete {
    println LikeletUtils.print_green("=================================================")
    println LikeletUtils.print_green("Cheers! RNAseq Pipeline from SYSUCC run Complete!")
    println LikeletUtils.print_green("=================================================")
    //email information
    if (params.mail) {
        recipient = params.mail
        def subject = 'My RNAseq-SYSUCC execution'

       def  msg = """\
            RNAseq-SYSUCC execution summary
            ---------------------------
            Your command line: ${workflow.commandLine}
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
        
            """.stripIndent()

        sendMail(to: recipient, subject: subject, body: msg)
    }
}
workflow.onError {

    println LikeletUtils.print_yellow("Oops... Pipeline execution stopped with the following message: ")+print_red(workflow.errorMessage)
}







def minimalInformationMessage() {

    println LikeletUtils.print_purple("You are running RNAseqPipe-SYSUCC with the following parameters:")
    println LikeletUtils.print_purple("Checking parameters ...")
    println LikeletUtils.print_yellow("=====================================")
    println LikeletUtils.print_yellow("Fastq file extension:           ") + print_green(params.read)
    println LikeletUtils.print_yellow("Single end :                    ") + print_green(params.singleEnd)
    println LikeletUtils.print_yellow("Strand specific condition:      ") + print_green(params.strand)
    println LikeletUtils.print_yellow("Output folder:                  ") + print_green(params.outdir)
    println LikeletUtils.print_yellow("STAR index path:                ") + print_green(params.star_index)
    println LikeletUtils.print_yellow("GTF path:                       ") + print_green(params.gene_gtf)
    println LikeletUtils.print_yellow("") + print_green(params.designfile)
    println LikeletUtils.print_yellow("Compare file path:              ") + print_green(params.comparefile)
    println LikeletUtils.print_yellow("=====================================")
    println "\n"
    // Minimal information message
    println LikeletUtils.print_green("-------------------------------------------------------------")
    println LikeletUtils.print_green("                       Checking Parameters                   ")
    println LikeletUtils.print_green("-------------------------------------------------------------")
    checkAnalysis("\tFastq file extension:           ",params.read)
    checkAnalysis("\tSingle end :                    ",params.singleEnd)
    checkAnalysis("\tStrand specific condition:      ",params.strand)
    checkAnalysis("\tOutput folder:                  ",params.outdir)
    checkAnalysis("\tSTAR index path:                ",params.star_index)
    checkAnalysis("\tGTF path:                       ",params.gene_gtf)
    checkAnalysis("\tDesign file  path:              ",params.designfile)
    checkAnalysis("\tCompare file path:              ",params.comparefile)
    println LikeletUtils.print_green("-------------------------------------------------------------")
}


def helpMessage(){

    println ''
    println LikeletUtils.print_purple('------------------------------------------------------------------------')
    println "RNAseqPipe_SYSUCC:  v$version"
    println LikeletUtils.print_purple('------------------------------------------------------------------------')
    println ''
    println LikeletUtils.print_yellow('Usage: ')
    println LikeletUtils.print_yellow('    The typical command for running the pipeline is as follows (we do not recommend users passing configuration parameters through command line, please modify the config.file instead):\n') 
            LikeletUtils.print_purple('       Nextflow run RNAseqPipe/main.nf ') +

            LikeletUtils.print_yellow('    General arguments:             Input and output setting') 
            print_parmeter('--inputdir <path> ','Path to input data(optional), current path default') 
            print_parmeter('--reads <*_fq.gz> ','Filename pattern for pairing raw reads, e.g: *_{1,2}.fastq.gz for paired reads') 
            print_parmeter('--outdir <path> ','The output directory where the results will be saved(optional), current path is default') 
            print_parmeter('Options: General options for run this pipeline') 
            print_parmeter('--designfile <file>' ,'A flat file stored the experimental design information ( required when perform differential expression analysis)') 
            print_parmeter('--comparefile <file>' ,'A flat file stored comparison information ( required when perform differential expression analysis, e.g )') 
            print_parmeter(' --singleEnd','Reads type, True for single ended ') 
            print_parmeter('--unstrand','RNA library construction strategy, specified for \'unstranded\' library ') 
            print_parmeter('--without_replicate ,'Specified when no replicates design involved, must provide  \'compare.txt\' file at the same time') 
            print_parmeter('--IDEA' ,'Run pre processing step for IDEA(http://idea.renlab.org) )')
            LikeletUtils.print_yellow('    References:                      If not specified in the configuration file or you wish to overwrite any of the references.') 
            print_parmeter('--fasta','Path to Fasta reference(required)') +
            print_parmeter('--gene_gtf' ,'An annotation file from GENCODE database in GTF format (required)\n') 
            LikeletUtils.print_yellow('    Other options:                   Specify the email and ') +
            print_parmeter('--mail' ,'email info for reporting status of your LncPipe execution  \n') +



    println '------------------------------------------------------------------------'
    println LikeletUtils.print_yellow('Contact information: zhaoqi@sysucc.org.cn')
    println LikeletUtils.print_yellow('Copyright (c) 2013-2018, Sun Yat-sen University Cancer Center.')
    println '------------------------------------------------------------------------'

}

def print_parameter(content, parameter){
    println LikeletUtils.print_cyan(LikeletUtils,addstringToalign(content, 30))+LikeletUtils.print_green(parameter)
}