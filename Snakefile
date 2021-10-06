configfile: "config.yaml"
import pandas as pd
import numpy as np

#ruleorder: trim > tobam > sort > feature_count> DGE

rule all:
      input:
          expand("galore/{SAMPLE}.r_1_val_1.fq.gz", SAMPLE=config['SAMPLE']), 
	  expand("galore/{SAMPLE}.r_2_val_2.fq.gz", SAMPLE=config['SAMPLE']), 
          expand("{SAMPLE}.bam", SAMPLE=config['SAMPLE']), 
          expand("{SAMPLE}.sorted.bam", SAMPLE=config['SAMPLE']), 
          expand("{treat}.counts.txt", treat= config['TREAT_NAME']),
          expand("{control}.counts.txt", control =config['CONTROL_NAME'])
rule trim: 
    input: 
       r1 = "{sample}.r_1.fq.gz",
       r2 = "{sample}.r_2.fq.gz"
    output: 
      "galore/{sample}.r_1_val_1.fq.gz",
      "galore/{sample}.r_2_val_2.fq.gz"
    conda: 'env/env-trim.yaml'
    shell: 
        """
         mkdir -p galore 
         mkdir -p fastqc 
         trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
        """ 

rule tobam:
    input:
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    params:
        genome=config['GENOME']
    output:
        "{sample}.bam"
    conda: 'env/env-align.yaml'
    shell:
        """
         bbmap.sh -Xmx24g in={input[0]} in2={input[1]} out={output} ref={params}
        """ 


rule sort:
    input: 
       "{sample}.bam"
    output:
       "{sample}.sorted.bam" 
    params: 
        "{sample}.tmp.sorted"
    log: 
        "{sample}.sorted.log" 
    conda: 'env/env-align.yaml'
    shell: 
       "samtools sort -T {params} -n -o {output} {input}"


rule feature_count:
     input: 
        expand("{sample}.sorted.bam", sample = config['TREAT']), 
        expand("{sample}.sorted.bam", sample = config['CONTROL'])
     params:
        config['STRAND'],
        config['GTF']
     output:
         expand("{treat}.counts.txt", treat =config['TREAT_NAME']),
         expand("{control}.counts.txt", control =config['CONTROL_NAME'])
     conda: 'env/env-feature.yaml'
     shell:
         """
          featureCounts -p -t exon -g gene_id -a {params[1]} -o {output[0]} {input[0]} -s {params[0]}
          featureCounts -p -t exon -g gene_id -a {params[1]} -o {output[1]} {input[1]} -s {params[0]}
         """ 
     
rule DGE: 
    input: 
       expand("{treat}.counts.txt", treat =config['TREAT_NAME']),
       expand("{control}.counts.txt", control =config['CONTROL_NAME'])
    output: 
         expand("{treat}.norm.txt", treat =config['TREAT_NAME']),
         expand("{control}.norm.txt", control =config['CONTROL_NAME'])
    conda: 'env/env-dge.yaml'
    shell: 
        "Rscript dge_genome.R {input[0]} {input[1]}" 
