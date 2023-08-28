#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//Modules
include { 
    BAMTOOLS_SPLIT_BAM;
    PICARD_MERGE_BAM;
} from '../modules/modules.nf'


//List of chromosomes(chr): from 1 to 23, X, Y, MT 
def chr = (1..23).collect{ 'chr' + it.toString() }
def xymt = ['chrX', 'chrY', 'chrMT']
chr.addAll(xymt)



// SPLIT BAM BY CHROMOSOME AND MERGE
workflow SPLIT_MERGE {
    take:
        bam_ch      // One bam as one element                            


    main:
        // Split bam by chromosome or scaffold (remainders are grouped as 'unmapped')
        bam_split_by_chr = BAMTOOLS_SPLIT_BAM(bam_ch) | flatten
       

        // Group bam by scaffold 
        bam_fList = bam_split_by_chr.map {                      //tuple structure: [ scaffold name, [ bam files ] ]                  
            file -> tuple(file.name.toString().tokenize('_').last().tokenize('.')[0..-2].join('.'), file) 
        }.groupTuple(by: 0)


        // Collect bam aligned to chromosomes
        bamList_chr = bam_fList.filter{ it[0] in chr }  
        

        // Merge bam by chromosome
        merged_bList = PICARD_MERGE_BAM(bamList_chr)       


    emit:
        merged_bList       // tuple structure: [ chromosome name, bam file ]  

}






