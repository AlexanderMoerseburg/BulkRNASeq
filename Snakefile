configfile: "config.yaml"

with open(config['TREAT']) as fp:
    TREAT = fp.read().splitlines()
with open(config['CONTROL']) as fp:
    CONTROL = fp.read().splitlines()

print(TREAT) 
print(CONTROL)

rule all:
      input:
          expand("{sample}.bam", sample = TREAT),
          expand("{sample}.bam", sample = CONTROL),
          expand("{sample}.sorted.bam", sample = TREAT),
          expand("{sample}.sorted.bam", sample = CONTROL),
          expand("{sample}.{group}.txt", sample = TREAT, group = config['TREAT_NAME']),
          expand("{sample}.{group}.txt", sample = CONTROL, group = config['CONTROL_NAME']),
          expand("{treat}_{control}_cpm.csv", treat =config['TREAT_NAME'],control =config['CONTROL_NAME']),
          expand("{treat}_{control}_dge.csv", treat =config['TREAT_NAME'],control =config['CONTROL_NAME'])

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
        "{sample}.sorted.bam"
     params:
        config['STRAND'],
        config['GTF'], 
        file = "tmp.txt"
     output:
         "{sample}.{group}.txt"
     conda: "env/env-feature.yaml"
     shell:
         """
          featureCounts -p -t exon -g gene_id -a {params[1]} -o {params[2]} {input[0]} -s {params[0]} 
          tail -n +3 {params[2]} | cut -f1,7 > {output[0]}
         """

rule DGE: 
    params: 
         treat = config['TREAT_NAME'], 
         control = config['CONTROL_NAME'], 
         N = config['N']
    output:
         expand("{treat}_{control}_cpm.csv", treat =config['TREAT_NAME'],control =config['CONTROL_NAME']),
         expand("{treat}_{control}_dge.csv", treat =config['TREAT_NAME'],control =config['CONTROL_NAME'])
    conda: 'env/env-dge.yaml'
    shell: 
         """
          Rscript scripts/dge_genome.R {params[0]} {params[1]} {params[2]}  
         """ 
