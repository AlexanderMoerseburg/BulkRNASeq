configfile: "config.yaml"
import pandas as pd

TREAT = pd.read_csv(config['TREAT'], header = 0)
CONTROL = pd.read_csv(config['CONTROL'], header = 0)

rule all:
      input:
          expand("{treat}.counts.txt", treat= config['TREAT_NAME']),
	  expand("{control}.counts.txt", control= config['CONTROL_NAME'])

if config['PAIRED']: 
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
else: 
     rule trim:
       input:
           "{sample}.fq.gz",
       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input} 
           """
if config['PAIRED']: 
    rule tobam:
         input:
              r1 = "galore/{sample}.r_1_val_1.fq.gz",
              r2 = "galore/{sample}.r_2_val_2.fq.gz"
         params:
              genome=config['GENOME'],
              mem = config['MEMORY'],
         output:
              "{sample}.bam"
         conda: 'env/env-align.yaml'
         shell:
             """
             bbmap.sh {params.mem} in={input[0]} in2={input[1]} out={output} ref={params.genome}
             """ 
else: 
      rule tobam:
         input:
              "galore/{sample}_trimmed.fq.gz",
         params:
              genome=config['GENOME'],
              mem = config['MEMORY'],
         output:
              "{sample}.bam"
         conda: 'env/env-align.yaml'
         shell:
             """
             bbmap.sh {params.mem} in={input} out={output} ref={params.genome}
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
       """samtools sort -T {params} -n -o {output} {input}"""

rule feature_count:
     input: 
        expand("{groupA}.sorted.bam" , groupA = TREAT), 
        expand("{groupB}.sorted.bam" , groupB = CONTROL)
     params:
        config['STRAND'],
        config['GTF']
     output:
         expand("{groupA}.counts.txt", groupA = config['TREAT_NAME']),
         expand("{groupB}.counts.txt", groupB = config['CONTROL_NAME'])
     conda: "env/env-feature.yaml"
     shell:
         """
          featureCounts -p -t exon -g gene_id -a {params[1]} -o {output[0]} {input[0]} -s {params[0]}
          featureCounts -p -t exon -g gene_id -a {params[1]} -o {output[1]} {input[1]} -s {params[0]}
        """

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
        "Rscript scripts/dge_genome.R {input[0]} {input[1]}" 
"""
