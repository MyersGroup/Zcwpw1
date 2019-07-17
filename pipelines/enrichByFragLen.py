

# To test if enrichment is due to specific lengths of fragments

#snakemake --cores 12 --snakefile pipelines/enrichByFragLen.py -npr

#snakemake --snakefile pipelines/enrichByFragLen.py --dag | dot -Tsvg > pipelines/enrichByFragLen_dag.svg

# requires:
# ../single-cell/sequencing/metadata/hg38_sizes.chrom
# MAPeakCaller/Fragment_Position_{sample}.sorted.bed
# deeptools/beds/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.bed

rule all:
  input:
    "byfraglength/WTCHG_538916_223180_vs_WTCHG_538916_221156_enrichByFragLen.pdf",
    "byfraglength/fraglen_WTCHG_538916_223180_byEnrichQ.pdf",
    "byfraglength/fraglen_WTCHG_538916_223180_221156_Quantised.pdf"


rule sortByFragLength:
  input:
    "MAPeakCaller/Fragment_Position_{sample}.sorted.bed"
  output:
    "byfraglength/{sample}.bed"
  threads:
    2
  shell:
    """
    awk -v OFS='\t' '{{print $0, $3-$2; }}' {input} | sort -n -k4 -S5G --parallel=2 > {output}
    """

rule quantiseBeds:
  input:
    "byfraglength/{sample}.bed"
  output:
    expand("byfraglength/{{sample}}_0{quantile}", quantile=[0, 1, 2, 3, 4])
  params:
    "byfraglength/{sample}_"
  shell:
    """
    split -d -n l/5 additional-suffix=.bed {input} {params}
    """

rule totalBases:
  input:
    "byfraglength/{sample}.bed"
  output:
    "byfraglength/{sample}.total"
  shell:
    "awk 'BEGIN{{size=0;}}{{size = size + $3 - $2;}}END{{print size;}}' {input} > {output}"


rule bedtograph:
  #
  # Create converage for IGV visualisation
  #
  input:
    bed="byfraglength/{sample}.bed",
    sizes="../single-cell/sequencing/metadata/hg38_sizes.chrom"
  output:
    "byfraglength/{sample}.bedgraph"
  threads:
    2
  shell:
    """
    LC_ALL=C sort {input.bed} -k1,1V -k2,2n -S5G --parallel=2 | 
    bedtools genomecov -bg -i stdin -g {input.sizes} | 
    LC_ALL=C sort -k1,1 -k2,2n -S5G --parallel=2 > {output}
    """


rule graphtobigwig:
  #
  # Compress bedgraph to Bigwig for IVG visualisation
  #
  input:
    graph="byfraglength/{sample}.bedgraph",
    sizes="../single-cell/sequencing/metadata/hg38_sizes.chrom"
  output:
    "byfraglength/{sample}.bigWig"
  threads:
    1
  shell:
    "bedGraphToBigWig {input.graph} {input.sizes} {output}"

rule convertToBED6:
  input:
    "deeptools/beds/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.bed"
  output:
    "byfraglength/hP9N_motifCenterStranded.bed"
  shell:
    """
    awk -v OFS='\t' '{{print $1,$2,$3,"0","0",$6; }}' {input} > {output}
    """

rule sumDepth:
  input:
    bigwig="byfraglength/{sample}.bigWig",
    regions="byfraglength/hP9N_motifCenterStranded.bed"
  output:
    "byfraglength/{sample}.profile"
  shell:
    "bwtool aggregate 2000:2000 {input.regions} {input.bigwig} {output} -fill=0"

rule plotAverages:
  input:
    sample_profile=expand("byfraglength/{{sample}}_0{quantile}.profile", quantile=[0, 1, 2, 3, 4]),
    input_profile=expand("byfraglength/{{input}}_0{quantile}.profile", quantile=[0, 1, 2, 3, 4]),
    totals=expand("byfraglength/{{sample}}_0{quantile}.total", quantile=[0, 1, 2, 3, 4])
  output:
    "byfraglength/{sample}_vs_{input}_enrichByFragLen.pdf"
  shell:
    "Rscript pipelines/plot_averages.R {output} {input.sample_profile} {input.input_profile}"


rule getOverlapingFragments:
  # get fragments that overlap different quantiles of ZcwCVC peaks
  input:
    fragment_position = "MAPeakCaller/Fragment_Position_{sample}.sorted.bed",
    peaks = "deeptools/beds/SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL.bed{quantile}"
  output:
    "byfraglength/Fragment_Position_{sample}_{quantile}.bed"
  shell:
    """
    bedtools intersect -a {input.fragment_position} -b {input.peaks} -wa > {output}
    """

rule plotFragLen:
# distribution of fragment lengths from each of those peak quantiles
  input:
    expand("byfraglength/Fragment_Position_{{sample}}_0{quantile}.bed", quantile=[0, 1, 2, 3]) # NB quartile not quintile!
  output:
    "byfraglength/fraglen_{sample}_byEnrichQ.pdf"
  shell:
    """
    Rscript pipelines/fragment_lengths.R {output} {input}
    """

rule plotQuantisedFragLen:
# distribution of fragment lengths from each of those peak quantiles
  input:
    sample=expand("byfraglength/{{group}}_{{sample}}_0{quantile}.bed", quantile=[0, 1, 2, 3, 4]),
    input=expand("byfraglength/{{group}}_{{input}}_0{quantile}.bed", quantile=[0, 1, 2, 3, 4])
  output:
    "byfraglength/fraglen_{group}_{sample}_{input}_Quantised.pdf"
  shell:
    """
    Rscript pipelines/fragment_lengths.R {output} {input.sample} {input.input}
    """

#Rscript fragment_lengths.R fraglen_WTCHG_538916_221156_Quantised.pdf byfraglength/WTCHG_538916_221156_0*.bed
#Rscript fragment_lengths.R fraglen_WTCHG_538916_Quantised.pdf byfraglength/WTCHG_538916_*_0*.bed

