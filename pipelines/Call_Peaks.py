# snakemake --cores 15 --snakefile Snakefile_peaks -npr




configfile: "pipelines/config.yml"
include: "sample_pairings.py"

from os.path import join
METADATA_DIR = config["metadata_dir"]
GENOME = 'hg38'



rule all:
  input:
    ["peaks/SingleBasePeaks.{chip}_vs_{control}.p0.000001.sep250.ALL.bed".format(
        chip = chip_id,
        control=control_id) for chip_id, control_id, name in sample_pairings],
    #["deeptools/bigwigs/{chip}_vs_{control}_CPM_log2.bw".format(
    #    chip = chip_id,
    #    control=control_id) for chip_id, control_id, name in sample_pairings]
    #["deeptools/bigwigs/SBPsample/{chip}_vs_{control}_CPM_log2.bw".format(
    #    chip = chip_id,
    #    control=control_id) for chip_id, control_id in sample_pairings]


def get_chips(wildcards):
  """
  When the input is required to be two file, return either the two biological replicates, or the two pseudoreplicates
  """
  if "_AND_" in wildcards.chip:
    two_samples = wildcards.chip.split("_AND_")
    return(list([f'FragPos/Fragment_Position_{two_samples[0]}.sorted.bed',
                f'FragPos/Fragment_Position_{two_samples[1]}.sorted.bed']))
  else:
    return(list([f'FragPos/Fragment_Position_{wildcards.chip}.sorted.bed.PR1',
                f'FragPos/Fragment_Position_{wildcards.chip}.sorted.bed.PR2']))

def split_chip_bam(wildcards):
  """
  When the input is required to be two file, split the wildcard on AND
  """
  if "_AND_" in wildcards.chip:
    two_samples = wildcards.chip.split("_AND_")
    return(list([f'filtered/{two_samples[0]}.bam',
                f'filtered/{two_samples[1]}.bam']))
  else:
    #print(wildcards.chip)
    return(list([f'NA',
                f'NA']))


#def get_1chip(wildcards):
#  """
#  When the input is required to be a single file, but in some cases there are two replicates, choose the first replicate
#  """
#  if "_AND_" in wildcards.chip:
#    two_samples = wildcards.chip.split("_AND_")
#    return(list([f'filtered/{two_samples[0]}.bam'])) # use first replicate
#  else: # no real replicates
#    return(list([f'filtered/{wildcards.chip}.bam']))

rule mergeFragPos:
  """
  Merge to fragment position bed files e.g. for replicates
  """
  input:
    a="FragPos/Fragment_Position_{sampleA}.sorted.bed",
    b="FragPos/Fragment_Position_{sampleB}.sorted.bed"
  output:
    sample = "FragPos/Fragment_Position_{sampleA}_AND_{sampleB}.sorted.bed",
  threads:
    5
  shell:
    """
    cat {input.a} {input.b} | LC_ALL=C sort -k1,1V -k2,2n -S5G --parallel=5 > {output}
    """

rule estconst:
  input:
    chips = get_chips ,
    control="FragPos/Fragment_Position_{control}.sorted.bed",
    chrsizes=join(METADATA_DIR, "hg38_sizes.chrom"),
  output:
    c="peaks/Constants.{chip}_vs_{control}.tsv",
  threads:
    3
  shell:
    """
    Rscript EstimateConstants.R \
	    peaks/ \
	    {input.chrsizes} \
	    {input.chips[0]} \
	    {input.chips[1]} \
      {input.control} \
      22 \
      {output.c}
    """


#rule call_peaks:
#  input:
#    chips = get_chips ,
#    control="FragPos/Fragment_Position_{control}.sorted.bed",
#    chrsizes=join(METADATA_DIR, "hg38_sizes.chrom")
#  output:
#    c="peaks/Constants.{chip}_vs_{control}.tsv",
#    p="peaks/SingleBasePeaks.{chip}_vs_{control}.p0.000001.sep250.ALL.bed"
#  threads:
#    15
#  shell:
#    """
#    sh FragPos.sh \
#	    --intermediates peaks/ \
#	    --chrsizes {input.chrsizes} \
#	    --outconstfile {output.c} \
#	    --outpeakfile {output.p} \
#	    -a {input.chips[0]} \
#	    -b {input.chips[1]} \
#     -i {input.control} \
#      --autosomes 22 \
#      --pthresh 0.000001 \
#      --peakminsep 250
#    """


rule call_peaks:
  input:
    chips = get_chips ,
    control="FragPos/Fragment_Position_{control}.sorted.bed",
    chrsizes=join(METADATA_DIR, "hg38_sizes.chrom"),
    consts="peaks/Constants.{chip}_vs_{control}.tsv"
  output:
    p="peaks/SingleBasePeaks.{chip}_vs_{control}.p0.000001.sep250.ALL.bed"
  threads:
    15
  shell:
    """
    Rscript DeNovoPeakCalling-SingleBase.R \
        peaks/ \
        {input.chrsizes} \
        {output.p} \
        {input.chips[0]} \
	      {input.chips[1]} \
        {input.control} \
        22 \
        0.000001 \
        250 \
        {input.consts}
    """


rule bamMerge:
  input:
    chips = split_chip_bam,
  output:
    "filtered/{chip}.bam"
  threads:
    1
  shell:
    """
    samtools merge {output} {input.chips[0]} {input.chips[1]}
    """

rule index_bam:
  input:
    "filtered/{sample}.bam"
  output:
    "filtered/{sample}.bam.bai"
  threads:
    5
  shell:
    "samtools index -@ 5 {input}"


rule bamCompare:
  input:
    chip="filtered/{chip}.bam",
    chip_indexes="filtered/{chip}.bam.bai",
    control="filtered/{control}.bam",
    control_indexes="filtered/{control}.bam.bai",
    black=join(METADATA_DIR, "{GENOME}.blacklist_notoverlaping.bed".format(GENOME=GENOME))
  output:
    "deeptools/bigwigs/{chip}_vs_{control}_CPM_log2.bw"
  threads:
    15
  shell:
    """
    bamCompare -b1 {input.chip} \
      -b2 {input.control} \
      --scaleFactorsMethod None \
      --normalizeUsing CPM \
      --operation log2 \
      --extendReads \
      --centerReads \
      -o {output} \
      --blackListFileName {input.black} \
      -p 15
    """

rule bamCompareSBP:
  input:
    chip="filtered/{chip}.bam",
    control="filtered/{control}.bam",
    black=join(METADATA_DIR, "{GENOME}.blacklist_notoverlaping.bed".format(GENOME=GENOME))
  output:
    "deeptools/bigwigs/SBPsample/{chip}_vs_{control}_CPM_log2.bw"
  threads:
    5
  shell:
    """
    bamCompare -b1 {input.chip} \
      -b2 {input.control} \
      --scaleFactorsMethod None \
      --normalizeUsing CPM \
      --operation log2 \
      --extendReads \
      --centerReads \
      --binSize 1 \
      -r chr10 \
      -o {output} \
      --blackListFileName {input.black} \
      -p 5
    """

rule deoverlap_blacklist:
  input:
    join(METADATA_DIR, "{GENOME}.blacklist.bed".format(GENOME=GENOME))
  output:
    join(METADATA_DIR, "{GENOME}.blacklist_notoverlaping.bed".format(GENOME=GENOME))
  shell:
    "sort -k1,1 -k2,2n {input} | bedtools merge -i stdin > {output}"


