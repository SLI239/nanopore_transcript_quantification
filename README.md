# Nanopore RNA-seq Analysis: Transcript Reconstruction and Quantification


## Introduction
**nanopore_transcript_quantification** is a bioinformatics best-practice analysis pipeline to reconstruct and quantify transcripts using Oxford Nanopore RNA-seq dataset.

## Requirements
  - Slurm
  - Miniconda3

## Tools
  - [Minimap2(v2.24)](https://github.com/lh3/minimap2)
  - [Samtools(v1.16.1)](http://www.htslib.org/doc/samtools.html)
  - [Bamtools(v2.5.2)](https://github.com/pezmaster31/bamtools)
  - [Gatk4(v4.3.0.0)](https://gatk.broadinstitute.org/hc/en-us/sections/9570256638747-4-3-0-0-Current)
  - [Stringtie(v2.2.1)](https://github.com/gpertea/stringtie)
  - [FeatureCounts(subread v2.0.3)](https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf)

## Input
- Gzipped fastq (.fastq.gz): Directory path needed
- Reference genome (.fa): File path needed
- Genome annotation (.gtf): File path needed 

## Output
Output is grouped by chromosome
- **chromosome name**_FCNTS_GENE.txt: read counts (FeatureCounts)
- **chromosome name**_FCNTS_GENE.txt.jcounts: Counts of reads supporting each exon-exon junction (FeatureCounts)
- **chromosome name**_FCNTS_GENE.txt.summary: summary of counting results (FeatureCounts)
- **chromosome name**.gtf: genome annotations (StringTie)

## Quick Start
1. Git clone this repository 
2. Activate Conda environment (If you don't have it already, create one and install packages as below)
```
conda install -c bioconda nextflow=22.04.5
conda install -c conda-forge mamba=1.3.1
``` 
3. Open 'nextflow.config' and edit 'params' section 
```
  1) fastq_dir: Path of gzipped fastq (.fastq.gz) directory
  2) out_dir: Path of directory to store output
  3) trace_dir: Path of directory to store pipeline execution information
  4) ref_genome: Path of reference genome file (.fa) 
  5) genome_annot: Path of genome annotation file (.gtf)
     - The version needs to match with the reference genome
  6) rna_type: Long-read RNA-seq technology, 'direct_rna' or 'cdna' (case-insensitive)
     - direct_rna: Direct sequencing of native RNA strands
     - cdna: Traditional RNA sequencing technology
  7) num_of_top_secondary_mapping: Minimap2 -N option
  8) min_chaining_score_ratio: Minimap2 -p option 
  9) fcnts_strand: FeatureCounts -s option 
  10) min_overlapping_bases: FeatureCounts --minOverlap option 
  11) feature_type: FeatureCounts -t option (feature type)
      - 'exon' or 'transcript', separated by ',' (no space) when multiple types are provided 
```
4. Run this pipeline  
```
nextflow run main.nf
```

## Note
 This pipeline was tested with RNA-seq generated using 1D sequencing technology 

## Credits
 **nanopore_transcript_quantification** was written by Songeun Lee (songeunlee2@gmail.com)
