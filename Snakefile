configfile: "config.yaml"


rule all:
      input:
          expand("galore/{SAMPLE}.r_1_val_1.fq.gz", SAMPLE=config['SAMPLE']), 
	  expand("galore/{SAMPLE}.r_2_val_2.fq.gz", SAMPLE=config['SAMPLE']), 
          expand("{SAMPLE}.sam", SAMPLE=config['SAMPLE']), 
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
    shell: 
        """
         mkdir -p galore 
         mkdir -p fastqc 
         trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
        """ 

rule tosam:
    input:
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    params:
        genome = config['GENOME']
    output:
        "{sample}.sam"
    shell:
        "bowtie2 -x {params.genome} -1 {input.r1} -2 {input.r2} -S {output}"


rule tobam:
      input:
          "{sample}.sam"
      output:
          "{sample}.bam"
      log:
          "{sample}.bam.log"
      shell:
          "samtools view {input[0]} -S -b > {output[0]}"

rule sort:
    input: 
       "{sample}.bam"
    output:
       "{sample}.sorted.bam" 
    params: 
        "{sample}.tmp.sorted"
    log: 
        "{sample}.sorted.log" 
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
     shell:
         """
          featureCounts -p -t exon -g gene_id -a {params[1]} -o {output[0]} {input[0]} -s {params[0]}
          featureCounts -p -t exon -g gene_id -a {params[1]} -o {output[1]} {input[1]} -s {params[0]}
         """ 

        
