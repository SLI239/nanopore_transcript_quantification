#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//Modules
include { 
    MINIMAP2_CREATE_MMI;
    MINIMAP2_ALIGN_READ;
    SAMTOOLS_SAM_TO_BAM;
} from '../modules/modules.nf'



// ALIGN READS AND GENERATE BAM 
workflow ALIGN {
    take:
        fastq_ch        // One fastq.gz as one element                            


    main:
        fq_fList = fastq_ch.map {                   // tuple structure: [ fastq name, fastq file ]
            file -> tuple(file.name.toString().tokenize('.')[0..-3].join('.'), file) 
        }
    

        // Referece genome indexing
        refIndex = MINIMAP2_CREATE_MMI(params.ref_genome)

        
        // Read alignment using reference index
        fq_refIndex = fq_fList.combine(refIndex)

        sam_ch = MINIMAP2_ALIGN_READ(fq_refIndex) 
        
        
        // Converting sam to bam 
        bam_ch = SAMTOOLS_SAM_TO_BAM(sam_ch)


    emit:
        bam_ch      // One bam as one element
        
}