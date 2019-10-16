
from os.path import join

configfile: '../pipelines/config.yml'

METADATA_DIR = config["metadata_dir"]

sample, strand, gz = glob_wildcards("ssDNA_{sample}_rep1_type1_filtered_only_rmdup.chrALL.{strand}prime.{gz}")

rule all:
    input:
        "DMC1_SSDS.pdf"

# Manual Step:
# combine seperate chromosome files into one
# cat ssDNA_ZCWPW1_HOM_260619_rep1_type1_filtered_only_rmdup.chr*.5prime.bedgraph.gz > ssDNA_ZCWPW1_HOM_260619_rep1_type1_filtered_only_rmdup.chrALL.5prime.bedgraph.gz
# 
# cat ssDNA_ZCWPW1_HOM_260619_rep1_type1_filtered_only_rmdup.chr*.3prime.bedgraph.gz > ssDNA_ZCWPW1_HOM_260619_rep1_type1_filtered_only_rmdup.chrALL.3prime.bedgraph.gz
# 
# cat WT/ssDNA_B6_Sample1_Brick2012_rep1_type1_filtered_only_rmdup.chr*.3prime.bedgraph > WT/ssDNA_B6_Sample1_Brick2012_rep1_type1_filtered_only_rmdup.chrALL.3prime.bedgraph 
# 
# cat WT/ssDNA_B6_Sample1_Brick2012_rep1_type1_filtered_only_rmdup.chr*.5prime.bedgraph > WT/ssDNA_B6_Sample1_Brick2012_rep1_type1_filtered_only_rmdup.chrALL.5prime.bedgraph 


rule gzip:
  input:
    "{sample}"
  output:
    "{sample}.gz"
  shell:
    "gzip -k {input}"


rule clean:
  # remove odd chromosomes
  # reomve Mitochondria, Sex
  # convert chr1 to 1
  input:
    "ssDNA_{sample}_rep1_type1_filtered_only_rmdup.chrALL.{strand}prime.bedgraph.gz"
  output:
    "ssDNA_{sample}_rep1_type1_filtered_only_rmdup.chrALLclean.{strand}prime.bedgraph"
  shell:
    """
    zgrep -v '_' {input} | grep -vP 'M|X|Y' | sed 's/chr//' > {output}
    """


rule bedtoBigWig:
  #
  # compress to bigwig for faster downstream processing
  #
  input:
    "ssDNA_{sample}_rep1_type1_filtered_only_rmdup.chrALLclean.{strand}prime.bedgraph"
  output:
    "{sample}_{strand}prime.bigWig"
  params:
    join("../",METADATA_DIR, "mm10_sizes.chrom")
  shell:
    """
    bedGraphToBigWig {input} {params} {output}
    """


rule averageProfile:
  input:
    sample="{sample}_{strand}prime.bigWig",
    b6="B6.bed",
    ko="KO.bed"
  output:
    b6="{sample}_{strand}prime_atB6.tsv",
    ko="{sample}_{strand}prime_atKO.tsv"
  params:
    a=5000,
    b=2000
  shell:
    """
    bwtool aggregate {params.a}:{params.b} {input.b6} {input.sample} {output.b6} -fill=0
    bwtool aggregate {params.a}:{params.b} {input.ko} {input.sample} {output.ko} -fill=0
    """


rule beds:
  # create bed file from composite file form Anjali
  #1 = chr (20 = X?)
  #4 = allele (B6, KO, or UNK)
  #5 = hshared
  #10 = motif center pos
  input:
    b6="B6_composite.txt",
    ko="KO.txt"
  output:
    ko="KO.bed",
    b6="B6.bed"
  shell:
    """
    awk -v OFS='\t' '$10=sprintf("%.0f",$10) {{if($1 != "20" && $5=="0" && $4=="B6") print $1,$10,$10,"0","0","+";}}'  {input.b6} > {output.b6}
    awk -v OFS='\t' '$10=sprintf("%.0f",$10) {{if($1 != "20" && $5=="0" && $4=="B6") print $0;}}'  B6_composite.txt > B6_composite_filtered.bed
    awk -v OFS='\t' '{{if($1 != "20") print $1,$2,$2,"0","0","+";}}'  {input.ko} | tail -n +2 > {output.ko}
    """

rule plotdmc1:
  input:
    expand("{sample}_{strand}prime_atB6.tsv", sample=set(sample), strand=set(strand)),
    expand("{sample}_{strand}prime_atKO.tsv", sample=set(sample), strand=set(strand))
  output:
    "DMC1_SSDS.pdf"
  shell:
    "Rscript ../pipelines/plotDMC1.R"
