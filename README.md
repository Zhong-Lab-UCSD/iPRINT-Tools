# iPRINT-Tools
Data decoding pipeline to process iPRINT-seq fastq data

## Overview
_in vivo_ protein-protein interaction sequencing (iPRINT-seq) efficiently maps cell-wide protein-protein interactions _in vivo_. iPRINT generally includes: **A. Labeling - B. Encoding - C. Decoding** steps.<br /> 
Here, we distribute iPRINT-Tools, a standardized data processing pipeline for the decoding step to identify protein-protein interactions from fastq files of iPRINT-seq experiments. The workflow of iPRINT is described as follows:<br />
<p align="center">
  <img src="https://github.com/Zhong-Lab-UCSD/iPRINT-Tools/blob/main/iPRINT_pipeline.png" width="66%">
</p>
The schematic diagram below describes the various stages of the iPRINT-Tools pipeline, including pre-processing of the raw reads, alignment to the transcriptome, identification of chimeric read pairs and identification of protein-protein interactions.<br />
<br />
<p align="center">
  <img src="https://github.com/Zhong-Lab-UCSD/iPRINT-Tools/blob/main/iPRINT_Decoding_Workflow.PNG" width="100%">
</p>
- At the pre-processing stage, with raw read pairs from the sequencing library as input, linker and adapter sequences are first removed. Low-quality and too short reads are then removed to get processed read pairs. <br /> 
- At the alignment stage, the pre-processed read pairs are mapped to the target transcriptome separately to get mapped read pairs. <br /> 
- At the next stage, we identify chimeric read pairs from the mapped read pairs. We select read pairs whose two ends’ primary alignments are mapped to different protein-coding genes and further check their mapping qualities. The read pairs passing the quality checks above are further deduplicated to be identified as chimeric read pairs. <br /> 
- At the stage of protein-protein interactions identification, for each chimeric read pair, we apply various statistical tests and cutoffs, including chi-square test, an odds ratio cutoff and a positive read count cutoff to finally identify protein-protein interactions (PPIs). <br /> 

## Decoding Workflow
1. Raw read pairs from the iPRINT-seq experiment are present in `.fastq` files.
2. Cutadpt is applied to remove 3' linker sequences and 5' adapter sequences from the read pairs. 
3. Fastp is then applied to remove low-quality reads whose mean quality is lower than Q20 and too short reads whose length is shorter than 20 bp.
3. The remaining read pairs are output as pre-processed read pairs in `.fastq` files.
4. The pre-processed read pairs are mapped to the transcriptome with BWA separately. ‘-a’ option is enabled to keep all found alignments using the default threshold of BWA. This is used in the later filtering of potential homologous read pairs. 
5. The mapped read pairs are output in `.csv` file with aligned genes and transcriptome alignment information.
6. The transcriptome alignment information of mapped read pairs is utilized to select read pairs whose two ends’ primary alignments are mapped to different protein-coding genes. The selected read pairs are further checked to see if both ends have over 50% of their read bases match the reference transcriptome based on the CIGAR string and if both ends have no shared lesser alignments. 
7. The read pairs passing the quality checks are then deduplicated based on the external coordinates of their primary alignments.
8. The deduplicated read pairs are identified as chimeric read pairs from the library. Their read ids and alignment information of the primary alignment are output in `chimericReadPairs.csv`. 
9. Chi-square test is applied to the chimeric read pairs. Benjamini-Hochberg adjustment is applied to correct all the p-values. Gene pairs with an adjusted p-value less than 0.05 (default) and with an odds ratio larger than 1 (default) are kept. Gene pairs with mapped chimeric read pair count in the library larger than 4 (default) times the average number of mapped chimeric read pairs per gene pair in the positive library are kept. 
10. The kept gene pairs are output as protein-protein interactions in `proteinProteinInteractions.csv`.


## Software Requirements
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) (1.18 or later)
- [fastp](https://github.com/OpenGene/fastp) (0.23.2 or later)
- [bwa](https://github.com/lh3/bwa) (0.7.17-r1188 or later)
- [samtools](http://www.htslib.org/) (1.13 or later)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (v2.30.0 or later)
- Python 3.6 or later, the following python libraries are required:<br />
    - sys
    - collections
    - cigar
    - glob
    - scipy
    - rpy2
    - datetime
## Additional files required
**BWA Index of the transcriptome to be aligned**<br />
You will need to download or build the bwa index of the target transcriptome (Human [RefSeq GRCh38 transcriptome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/) or Mouse [RefSeq GRCm39 transcriptome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/)) for iPRINT-Tools to use. Using:<br />
`bwa index /path/to/your/transcriptome.fa`

**Transcript, gene and gene type dictionary file**<br />
You will also need a dictionary file that contains the information of transcript ids to their corresponding gene names/gene ids and corresponding gene types in a csv format with the first column being transcript ids, the second column being gene names and the third column being gene types. Here we provide an example dictionary file for Human [RefSeq GRCh38 genome](https://github.com/Zhong-Lab-UCSD/iPRINT-Tools/blob/main/refSeq_tx_gene_type_human.csv) or Mouse [RefSeq GRCm39 genome](https://github.com/Zhong-Lab-UCSD/iPRINT-Tools/blob/main/refSeq_tx_gene_type_mouse.csv)

## Usage
**Installation**
1. Clone the current github repository to your local machine. For example<br />
`git clone https://github.com/Zhong-Lab-UCSD/iPRINT-Tools`
2. Add the following path of the cloned directory to your `.bashrc` file<br />
`export PATH=$PATH:/home/path/to/iPRINT-Tools/bin`

**To execute iPRINT-Tools, run**
<pre><code>
iPRINT-Tools -a /path/to/read1.fastq
             -b /path/to/read2.fastq
             -i /path/to/bwaIndex/transcriptome.fa
             -o /path/to/outputDir
             -g /path/to/refSeq_tx_gene_type.csv
           
</code></pre>


**Required parameters**
<pre><code>
-a     |String, Path to read1 fastq file, fastq.gz also supported
-b     |String, Path to read2 fastq file, fastq.gz also supported
-o     |String, Path to output directory
-i     |String, Path to the bwa index of the target transcriptome
-g     |String, Path to transcript, gene and gene type dictionary file
    
</code></pre>

**Other parameters**
<pre><code>
-d     |Float, odds ratio cutoff used to identify protein-protein interactions, default=1
-p     |Float, false discovery rate cutoff used to identify protein-protein interactions, default=0.05
-c     |Float, read count cutoff coefficient used to identify protein-protein interactions, default=2
-j     |String, Job ID to be prepended to the output files and directories, optional, default=iPRINT-Tools"
-t     |Int, Number of working threads, default=2
-r     |Char, (T or F), removal of intermediate files or not, default=T
-h     |Print usage message"     
    
</code></pre>

## iPRINT-Tools Output
A variety of output files is created for each sample as it is run through the pipeline. The highest level of the output directory contains the following files and subdirectories:
<pre><code>
proteinProteinInteractions.csv    |a file that contains the identified protein-protein interactions from the sample
chimericReadPairs.csv             |a file that contains the read ids of the identified chimeric read pairs from the sample
summary.csv                       |a file that contains the summary statistics of running the sample with iPRINT-Tools
errorLog.txt                      |a file that contains error messages from the pipeline if any
processedFastq/                   |a directory that contains the pre-processed fastq files from the sample
alignment/mappedReadPairs.csv     |a file that contains the alignment information of all mapped read pairs from the sample
alignment/*/                      |subdirectories that contain the alignment files of the pre-processed fastq files from the sample
intermediateFiles/                |an optional directory contains all the intermediate files generated from running the pipeline, this directory only exists if the '-r' option is set to 'F' 
    
</code></pre>
