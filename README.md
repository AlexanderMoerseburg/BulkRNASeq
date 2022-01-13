Snakemake Workflow for RNASeq Analysis 
==========================================================================


This is an RNASeq snakemake pipeline written by Sherine Awad. 
The pipeline uses BBmap, featurecounts, and egdeR to do analysis at genome level. It also uses Salmon, tximport, and egdeR to do analysis at transcripts level.  
More tools and more options are being added to the pipeline. 

#### Edit config file 

The config file has several sections. You will need to edit the general section as follows: 


Edit the config file, choose if you have paired-end or single reads. For ppaired end read, change PAIRED entry in config file to TRUE (default). For single end reads, change PAIRED entry in config file to FALSE. 

Edit the treatment.tsv and control.tsv to have your treatments and control names.  You can edit TREAT and CONTROL entries in config file to change these file names. 

Edit the TREAT_NAME  and CONTROL_NAME in the config file to change treatment and control names. Default is stress and wildtype respecively. 

Edit MEMORY entry in config file for how much GB needed. 

Edit LEVEL entry in config file to GENOME or TRANSCRIPTS whether to run at genome level or transcripts level respectively. 
 

If you choose to run at the genome level, you will need to edit the genome level section with the GENOME, GTF, and INDEX for genome, gtf file and index location respetively. You can ignore these entries if you choose to run at transcripts level.  

If you choose to run at the transcripts level, you will need to edit the transcripts level section with TRANSCRIPTS, SALMON_INDEX, SALMON_LIBRARY, TX2GENE, and TXGTF for transcript, SALMON index location, library type as used in SALMON (SR, SF, etc), tx2gene file name and transcripts GTF. 

Currently, we are using Salmon only for transcripts level, more tools will be used in the near future. 

Refere to our Makefile to prepare Salmon Index. See Prepare Reference section below. 


#### Prepare Reference

Whether you will run at genome level or transcripts level, you can use our  Makefile to pull the reference if needed (See below). 

For the sake of reproducibility, we provide a makefile to pull the reference; HUMAN or Mouse. More references will be included.


For Human Ensembl GRCh37:


    make Homo_sapiens_Ensembl_GRCh37.tar.gz

For Mouse Ensembl GRCm38:


    make Mus_musculus_Ensembl_GRCm38.tar.gz

If you choose to run at transcripts level using Salmon, you can prepare Salmon index as below:: 

   
   make SalmonGRCh37Index

or 

  make SalmonGRCh38Index

Depending on which transcript you choosed. 

 

#### Run Snakemake pipeline 

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


 
