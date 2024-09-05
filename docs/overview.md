# :dog::dog::dog::notes: Fluffy - Hematology WGS and WTS pipeline :snake:
The Fluffy pipeline is designed to be able to process both DNA and RNA whole transcriptome/genome sequencing and therefore have two parallel tracks. The DNA track you can run either with a match normal or as tumor only sample. The RNA part is designed to primarily identify fusions. 
Here follows a brief overview of the pipelines:

## DNA
The DNA-part of the pipeline allows for both tumor only analysis and tumor with a match normal analysis. It generates SNV-calls and summarizes them into different panels and formats. It also does CNV calling using a combination of both CNVkit and GATK to then summarize it into one html-report. To identify SVs Manta is used. It also runs several different QC-programs which are then summarized in an MultiQC-report. 



![dag plot](includes/images/dna.svg){: style="height:100%;width:100%"}

## RNA
Fluffy can also process whole transcriptome data to identify potential fusions. It uses three different fusion callers 

![dag plot](includes/images/rna.svg){: style="height:100%;width:100%"}