#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//Modules
include { 
    STRINGTIE_CREATE_GTF;
    FCNTS_QUANT_TRANSCRIPT;
 } from '../modules/modules.nf'



// RECONSTRUCT AND QUANTIFY TRANSCRIPTS
workflow QUANT {
    take:
        merged_bList       // tuple structure: [ chromosome name, bam file ]                                                                        


    main:
        // Assemble alignments into transcripts using StringTie
        strT_fList = STRINGTIE_CREATE_GTF(merged_bList)


        // Quantify transcripts using FeatureCounts
        bam_gtf = merged_bList.combine(strT_fList, by:0)
        
        fcnts_fList = FCNTS_QUANT_TRANSCRIPT(bam_gtf)
    
    
    emit:
        fcnts_fList      // FeatureCounts files for each chromosome (tuple) 

}