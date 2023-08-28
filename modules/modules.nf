// Generate reference genome index 
process MINIMAP2_CREATE_MMI {

    input:
        val(ref_fasta)
    output:
        file('*.mmi')

    script: 
    """
    minimap2 -t $task.cpus -d ref.mmi $ref_fasta    

    """   
}


// Align reads using reference genome index
// Retain up to -N top secondary mappings, 
// if their chaining scores are higher than -p of their corresponding primary mappings
process MINIMAP2_ALIGN_READ {

    input:
        tuple val(fastq_name), file(fastq_reads), file(ref_index)
    output:
        tuple val(fastq_name), file('*.sam')

    script:   
    if( params.rna_type.toLowerCase() == 'direct_rna' )   
    """
    minimap2 -t $task.cpus \
    -ax splice -uf -k14 \
    -p $params.min_chaining_score_ratio \
    -N $params.num_of_top_secondary_mapping \
    $ref_index $fastq_reads > ${fastq_name}.sam 
        
    """
    else if( params.rna_type.toLowerCase() == 'cdna')
    """
    minimap2 -t $task.cpus \
    -ax splice:hq -uf \
    -p $params.min_chaining_score_ratio \
    -N $params.num_of_top_secondary_mapping \
    $ref_index $fastq_reads > ${fastq_name}.sam 

    """
    else
        throw new IllegalArgumentException("ERROR:UNKNOWN RNA-SEQ DATA TYPE($params.rna_type)\nCHOOSE ONE OF THE FOLLOWING OPTIONS: direct_rna or cdna")

}


// Generate bam
process SAMTOOLS_SAM_TO_BAM {

    input:
        tuple val(sam_name), file(sam)
    output:
        file('*.bam')

    script: 
    """
    samtools view -bh -@ $task.cpus ${sam} > ${sam_name}.bam 

    """
}


// For each bam, split reads by chromosome 
process BAMTOOLS_SPLIT_BAM {

    input:
        file(bam)
    output:
        file('*')

    script: 
    """  
    bamtools split -in $bam -reference 

    """   
}


// Merge bam by chromosome using picard
process PICARD_MERGE_BAM {
    
    input:
        tuple val(chr), file(bamList_chr)
    output:
        tuple val(chr), file('*.bam')

    script:
    """
    input_bam=\$(ls *.bam | awk '{printf(" -I ", \$1); system("readlink -f " \$1)}')

    gatk --java-options '-Xmx${task.memory.toGiga()}g' \
    MergeSamFiles \
    \$input_bam \
    -O ${chr}.bam \
    --ASSUME_SORTED false \
    --USE_THREADING true \
    --SORT_ORDER coordinate \
    --CREATE_INDEX false
 
    """
}


// Create GTF by chromosome
process STRINGTIE_CREATE_GTF {
    publishDir "${params.out_dir}/${chr}", overwrite: false

    input:
        tuple val(chr), file(bam)
    output:
        tuple val(chr), file('*.gtf')

    script:
    """
    stringtie -L \
    -G $params.genome_annot \
    -o ${chr}.gtf \
    -v -p $task.cpus $bam
    
    """
}


// Count reads mapped to genomic features
process FCNTS_QUANT_TRANSCRIPT {
    publishDir "${params.out_dir}/${chr}", overwrite: false
     
    input:
        tuple val(chr), file(bam), file(gtf)                                            
    output:
        file('*')

    script:
    """
    featureCounts -T $task.cpus \
    -t $params.feature_type \
    -g gene_id \
    -a $gtf -F GTF \
    -G $params.ref_genome \
    -o ${chr}_FCNTS_GENE.txt \
    -C -J -L \
    --ignoreDup \
    --minOverlap $params.min_overlapping_bases \
    -s $params.fcnts_strand $bam 
    
    """
}

