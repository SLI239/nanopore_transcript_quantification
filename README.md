# Nanopore RNA-seq Analysis: Transcript Reconstruction and Quantification


## Introduction
**nanopore_transcript_quantification** is a bioinformatics best-practice analysis pipeline to reconstruct and quantify transcripts using Oxford Nanopore RNA-seq dataset.

It includes following tools. For detailed information, check out the links:
  - **minimap2 (v2.24)**: https://github.com/lh3/minimap2
  - **samtools (v1.16.1)**: http://www.htslib.org/doc/samtools.html
  - **bamtools (v2.5.2)**: https://github.com/pezmaster31/bamtools
  - **gatk4 (v4.3.0.0)**: https://gatk.broadinstitute.org/hc/en-us/sections/9570256638747-4-3-0-0-Current
  - **stringtie (v2.2.1)**: https://github.com/gpertea/stringtie
  - **featureCounts (subread v2.0.3)**: https://subread.sourceforge.net/ https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf


## Requirements
  - Slurm
  - Miniconda3


## Quick Start
1. Git clone this repository 
2. Activate Conda environment (If you don't have one already, create one and install packages as below)
```
conda install -c bioconda nextflow=22.04.5
conda install -c conda-forge mamba=1.3.1
``` 
3. Open nextflow.config and edit 'params' section 
```
  1) fastq_dir: Path of gzipped fastq (.fastq.gz) directory
  2) ref_genome: Path of reference genome file (.fa) 
  3) genome_annot: Path of genome annotation file (.gtf)
     - The version needs to match with the reference genome
  4) out_dir: Path of directory to store output
  5) trace_dir: Path of directory to store pipeline execution information
  6) rna_type: Long-read RNA-seq technology, 'direct_rna' or 'cdna' (case-insensitive)
     - direct_rna: Direct sequencing of native RNA strands
     - cdna: Traditional RNA sequencing technology
  7) num_of_top_secondary_mapping: Minimap2 -N option (number of top secondary mappings)
  8) min_chaining_score_ratio: 
     Minimap2 -p option (minimum ratio of secondary mapping chaining score to primary)  
  9) fcnts_strand: FeatureCounts -s option (strand-specific read counting)
     - 0 (unstranded; default), 1 (stranded) or 2 (reversely stranded)
  10) min_overlapping_bases: 
      FeatureCounts --minOverlap option (minimum number of overlapping bases in a read) 
      - default: 1
  11) feature_type: FeatureCounts -t option (feature type)
      - 'exon' or 'transcript', separated by ',' (no space) when multiple types are provided 
```
4. Run this pipeline  
```
nextflow run main.nf
```


## Input
- Gzipped fastq (.fastq.gz): Directory path required
- Reference genome (.fa): File path required
- Genome annotation (.gtf): File path required 


## Output
Output data is groupped by chromosome
- **chromosome name**_FCNTS_GENE.txt: read counts (FeatureCounts)
- **chromosome name**_FCNTS_GENE.txt.jcounts: Counts of reads supporting each exon-exon junction (FeatureCounts)
- **chromosome name**_FCNTS_GENE.txt.summary: summary of counting results (FeatureCounts)
- **chromosome name**.gtf: genome annotations (StringTie)


## Note
 This pipeline was tested with RNA-seq generated using 1D sequencing technology 


## Credits
 **nanopore_transcript_quantification** was written by Songeun Lee (songeunlee2@gmail.com)
