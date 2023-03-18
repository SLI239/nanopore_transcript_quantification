# Nanopore RNA-seq downstream analysis: Transcript reconstruction and quantification


## Introduction
**nanopore_transcript_quantification** is a bioinformatics best-practice analysis pipeline to reconstruct and quantify transcripts using Oxford Nanopore RNA-seq dataset.

It includes following tools. For detailed information, check out the links:
  - minimap2 (v2.24): https://github.com/lh3/minimap2
  - samtools (v1.16.1): http://www.htslib.org/doc/samtools.html
  - bamtools (v2.5.2): https://github.com/pezmaster31/bamtools
  - gatk4 (v4.3.0.0): https://gatk.broadinstitute.org/hc/en-us/sections/9570256638747-4-3-0-0-Current
  - stringtie (v2.2.1): https://github.com/gpertea/stringtie
  - featureCounts (subread v2.0.3): https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf


## Requirements
  - Slurm
  - Miniconda3


## Quick Start
1. git clone this repository 
2. Set up conda environment: 
    Activate a Conda environment (and install following packages if not installed already)
```
    conda install -c bioconda nextflow=22.04.5
    conda install -c conda-forge mamba=1.3.1
``` 
3. Open 'nextflow.config' and edit 'params' section 
```
    - fastq_dir: Path of fastq.gz directory
    - ref_genome: Path of reference genome file (.fa) 
    - genome_annot: Path of genome annotation file (.gtf); THE VERSION OF GENOME ANNOTATION NEEDS TO MATCH WITH THE REFERENCE GENOME  
    - out_dir: Path of directory to store output
    - trace_dir: Path of directory to store pipeline execution information
    - rna_type: Long-read RNA-seq technologies (case-insensitive); either "direct_rna" or "cdna"
                "direct_rna" refers to direct RNA and "cdna" refers to traditional full-length cDNA 
                For more information, check out https://github.com/lh3/minimap2
    - num_of_top_secondary_mapping: Minimap2 -N option; Check out https://github.com/lh3/minimap2
    - secondary_mapping_min_chaining_ratio: Minimap2 -p option; Check out https://github.com/lh3/minimap2
    - fcnts_strand: FeatureCounts -s option
                    Indicate if strand-specific read counting should be performed; 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
    - min_overlapping_bases: FeatureCounts --minOverlap option
                        Minimum number of overlapping bases in a read that is required for read assignment (default: 1) 
    - feature_type: FeatureCounts -t option
                    The feature type(s); 'exon' and/or 'transcript'
                    If more than one feature type is provided, they should be separated by ',' (e.g., "exon,transcript")  
```
4. Run this pipeline  
```
nextflow run main.nf
```


## Input
- Gzipped fastq: Directory path required
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
 nanopore_transcript_quantification was written by Songeun Lee (songeunlee2@gmail.com)
