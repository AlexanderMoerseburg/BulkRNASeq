Snakemake Workflow for RNASeq Analysis 
==========================================================================


This is an RNASeq snakemake pipeline written by Sherine Awad. 
The pipeline uses edgeR and BBmap. 
More tools and more options are being added to the pipeline. 

To run the pipeline, edit the config file to match your sample names, your reference genome then: 


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


## TODO 


#### Dry run 

for a dry run use:

    snakemake -j1 -n

and you can see the command printed on a dry run using:

    snakemake -j1 -n -p

#### Keep going option 


You can try the following to keep going if any issues happen, like no variants is found by one tool:

    snakemake -j1 --keep-going


## TODO 

1. Add DESEQ2 options 
2. Add Star and other aligners along with BBmap for more options 
3. Add mapping to transcriptomes like Salmon, etc. 


