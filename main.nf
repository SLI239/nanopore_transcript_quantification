#!/usr/bin/env nextflow
nextflow.enable.dsl=2



//Log header
log.info """\

              TRANSCRIPT RECONSTRUCTION & QUANTIFICATION

     - TOOLS: STRINGTIE (v2.2.1) AND FEATURECOUNTS (SUBREAD v2.0.3)
     - DATASET: OXFORD NANOPORE 1D RNA-SEQ (DIRECT RNA OR CDNA)
        
        
     * SETTINGS
     =====================================================================
     OPTION (RNA-SEQ TECHNOLOGY)      :       $params.rna_type
     INPUT (FASTQ.GZ DIRECTORY)       :       $params.fastq_dir 
     INPUT (REFERENCE GENOME FILE)    :       $params.ref_genome
     INPUT (GENOME ANNOTATION FILE)   :       $params.genome_annot
     OUTPUT (OUTPUT DIRECTORY)        :       $params.out_dir
     OUTPUT (RUN INFO DIRECTORY)      :       $params.trace_dir
        
"""



//Subworkflows
include { ALIGN } from './steps/align_fq_reads_step.nf'
include { SPLIT_MERGE } from './steps/split_merge_bam_step.nf'
include { QUANT } from './steps/cnt_transcript_step.nf'



//Workflow: main
workflow {     
    // Initial channel  
    fastq_ch = Channel.fromPath(params.fastq_dir + "/*.fastq.gz", checkIfExists: true)


    // 1) Align reads to the reference genome  
    ALIGN(fastq_ch)


    // 2) Merge mapped reads by chromsome
    SPLIT_MERGE(ALIGN.out)
    

    // 3) Transcript reconstruction and quantification 
    QUANT(SPLIT_MERGE.out)

}



//Print out workflow completion status
workflow.onComplete {
    log.info ( workflow.success ? "\nDone! \n" : "Oops.. something went wrong" )
}

 
