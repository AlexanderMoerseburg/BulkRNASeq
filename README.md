[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Snakemake Workflow for Bulk RNASeq  
===============================================


This is an BulkRNASeq snakemake pipeline written by Sherine Awad. 
The pipeline uses BBmap, featurecounts, and egdeR to do the analysis at genome level. It also uses Salmon, tximport, and egdeR to do the analysis at transcripts level.  
The pipeline works for any genome or transcripts. More tools and more options are being added to the pipeline. 

## Edit config file 

The config file has several sections. You will need to edit the general section as follows: 

Edit the general parameters: 

| Config Variable  | Description                      |
| -----------------| ---------------------------------|
| TREAT            | The file name that contains your treatment sample names. Default: treatment.tsv |
| CONTROL          | The file name that contains your control sample names. Default: control.tsv  |
| TREAT_NAME       | The name you wish to give to your treatment in figures. Default: stress|
| CONTROL_NAME     | The name you wish to give to your control in figures. Default: wildtype |
| WORKDIR          | your working directory |
| STRAND           | Strand: 1 or 2 according to feature counts definition of strands |
| N                | Number of replcates  |
| PAIRED           | True is your samples are paired-end, false if your samples are single-ended |
| Memory           | Memory available: "-Xmx40g" is the default for 40GB |
| LEVEL            | Set to  GENOME or TRANSCRIPT: Would you like to align to the genome or map to transcipts. Genome level vs transcripts level analysis |


Edit the treatment.tsv and control.tsv to have your treatments and control names.  You can edit TREAT and CONTROL entries in config file to change these file names.

The pipeline takes samples with a suffix 'r_1.fq.gz' and 'r_2.fq.gz' if the samples are paired. Or it takes samples with suffix 'fq.gz' if the samples is single-end reads.
Regardless your samples are paired or single-ended, samples names should be listed in treatment.tsv and control.tsv or their equivalent without the suffix.

Edit the TREAT_NAME  and CONTROL_NAME in the config file to change treatment and control names. Default is stress and wildtype respecively.


If you choose to run at genome level, edit the genome level parameters: 


| Config Variable  | Description                      |
| -----------------| ---------------------------------|
| GENOME           | Path to your genome.fa |
| GTF              | Path to your annotation file: genes.gtf   |
| INDEX            | Path to your index|


If you choose to run at the transcripts level, edit the transcripts level parameters: 

| Config Variable  | Description                      |
| -----------------| ---------------------------------|
| TRANSCRIPTS      | Path to your transcript.fa |
| SALMON_INDEX     | Path to your Salmon index  |
| SALMON_LIBRARY   | Library strand according to Salmon |
| TX2GENE          | Name of tx2gene.txt file for tximport |
| TXGTF            | Path to annotation file gtf | 


The tx2gene will generated automatically and saved into name given in TX2GENE entry in config file using our script txi2gene.R. 

Currently, we are using Salmon only for transcripts level, more tools will be used in the near future. 

You can used our Makefile to pull the reference genome and for Salmon Index preparation. See the section below: i

### Prepare Reference

For the sake of reproducibility, we provide a makefile to pull the reference; HUMAN or Mouse. But you can include any reference you wish to use. 
We are giving examples using Human and Mouse genome. 


For Human Ensembl GRCh37:


    make Homo_sapiens_Ensembl_GRCh37.tar.gz

For Mouse Ensembl GRCm38:


    make Mus_musculus_Ensembl_GRCm38.tar.gz

If you choose to run at transcripts level using Salmon, you need to pull the transcripts and GTF first: 


    make Homo_sapiens.GRCh38.cdna.all.fa.gz

Then::


    make Homo_sapiens.GRCh38.104.gtf.gz

Or use any other transcripts/GTF.  


Then you need to prepare Salmon index as below:: 

   
    make SalmonGRCh37Index

or: 

    make SalmonGRCh38Index


Depending on which transcripts you choosed. 


### Run Snakemake pipeline 

Once you edit the config file to match your needs, then:  


    snakemake -jn 

where n is the number of cores for example for 10 cores use:


    snakemake -j10 


For a dry run use: 
  
  
    snakemake -j1 -n 


and to print command in dry run use: 

  
    snakemake -j1 -n -p 

  
#### Use Conda 

For less frooodiness, to pull automatically the same versions of dependencies use:

    snakemake -jn --use-conda

This will pull the same versions of tools we used. Conda has to be installed in your system.

For example, for 10 cores:

    snakemake -j10 --use-conda


#### Dry run 

for a dry run use:

    snakemake -j1 -n

and you can see the command printed on a dry run using:

    snakemake -j1 -n -p


#### Keep going option 


You can try the following to keep going if any issues happen, like no variants is found by one tool:

    snakemake -j1 --keep-going


#### Collect some stats 

    snakemake -j 10 --keep-going --stats run.stats


### References 
1. Bushnell, B. (2014). BBMap: a fast, accurate, splice-aware aligner (No. LBNL-7065E). Lawrence Berkeley National Lab.(LBNL), Berkeley, CA (United States).
2. Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature methods, 14(4), 417-419.
3. Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140.
4. Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.
5. Luo, W., Friedman, M. S., Shedden, K., Hankenson, K. D., & Woolf, P. J. (2009). GAGE: generally applicable gene set enrichment for pathway analysis. BMC bioinformatics, 10(1), 1-17.
6. Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4.
7. Love, M. I., Soneson, C., & Robinson, M. D. (2017). Importing transcript abundance datasets with tximport. Dim Txi. Inf. Rep. Sample1, 1, 5.
 


