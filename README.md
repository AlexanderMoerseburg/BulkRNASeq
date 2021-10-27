Snakemake Workflow for RNASeq Analysis 
==========================================================================


This is an RNASeq snakemake pipeline written by Sherine Awad. 
The pipeline uses BBmap, featurecounts, and egdeR to do analysis at genome level. It also uses Salmon, tximport, and egdeR to do analysis at transcripts level.  
More tools and more options are being added to the pipeline. 

#### Edit config file 

The config file has several sections. You will need to edit the general section as follows: 


Edit the config file, choose if you have paired-end or single reads. Edit the treatment.tsv and control.tsv to have your treatments and control names. 

Then you you need to change the level whether to run at genome level or transcripts level. 

If you choose to run at the genome level, you will need to edit the genome level section with the reference, index, strand, etc. 

If you choose to run at the transcripts level, you will need to edit the transcripts level section with reference, index, strand, etc.

### Pull Reference

Whether you will run at genome level or transcripts level, you can use our  Makefile to pull the reference if needed (See below). 

For the sake of reproducibility, we provide a makefile to pull the reference; HUMAN or Mouse. More references will be included.


For Human Ensembl GRCh37:


    make Homo_sapiens_Ensembl_GRCh37.tar.gz

For Mouse Ensembl GRCm38:


    make Mus_musculus_Ensembl_GRCm38.tar.gz

####Run Snakemake pipeline 

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


 
