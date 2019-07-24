# snakemake --cores 15 --snakefile pipelines/Force_Call_Peaks.py -npr

configfile: "pipelines/config.yml"
include: "sample_pairings.py"


from os.path import join
METADATA_DIR = config["metadata_dir"]
GENOME = 'hg38'

# override sample_pairings to compute
# sample_pairings = [("SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139", "SRA_Altemose2015_SRR5627140"),# YFP Nterm
#                 ("SRA_Altemose2015_SRR5627141", "SRA_Altemose2015_SRR5627140")] # YFP Nterm H3K4me3

forcecallAT = config["forcecallAT"]

rule all:
  input:
    ["peaks/ForceCalledPeaks_{chip}_vs_{control}_AT_{at}.bed".format(
        chip = chip_id,
        control = control_id, at=forcecallAT) for chip_id, control_id, name in sample_pairings]


def get_chips(wildcards):
  if "_AND_" in wildcards.chip:
    two_samples = wildcards.chip.split("_AND_")
    return(list([f'MAPeakCaller/Fragment_Position_{two_samples[0]}.sorted.bed',
                f'MAPeakCaller/Fragment_Position_{two_samples[1]}.sorted.bed']))
  else:
    return(list([f'MAPeakCaller/Fragment_Position_{wildcards.chip}.sorted.bed.PR1',
                f'MAPeakCaller/Fragment_Position_{wildcards.chip}.sorted.bed.PR2']))


rule estconst:
  input:
    chips = get_chips ,
    control="MAPeakCaller/Fragment_Position_{control}.sorted.bed",
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



rule extend_peaks:
  input:
    forcepos="peaks/SingleBasePeaks.{forcecallAT}.p0.000001.sep250.ALL.bed",
    chrsizes=join(METADATA_DIR, "hg38_sizes.chrom")
  output:
    "peaks/SingleBasePeaks.{forcecallAT}.flank{flank}.p0.000001.sep250.ALL.bed"
  threads:
    1
  shell:
    """
    # bedtools doesn't like the header so use tail -n +2
    tail {input.forcepos} -n +2 | bedtools slop -i - -g {input.chrsizes} -b {wildcards.flank} > {output}
    """

rule call_peaks:
  input:
    chips = get_chips ,
    control="MAPeakCaller/Fragment_Position_{control}.sorted.bed",
    chrsizes=join(METADATA_DIR, "hg38_sizes.chrom"),
    forcepos= lambda wc: "forcepeaks/{forcecallAT}" if(".bed" in wc.forcecallAT) else "peaks/SingleBasePeaks.{forcecallAT}.p0.000001.sep250.ALL.bed",
    consts="peaks/Constants.{chip}_vs_{control}.tsv"
  output:
    p="peaks/ForceCalledPeaks_{chip}_vs_{control}_AT_{forcecallAT}.bed"
  threads:
    3
  shell:
    """
    Rscript ForceCallPeaks.R \
	    peaks/ \
	    {input.chips[0]} \
	    {input.chips[1]} \
      {input.control} \
      {input.consts} \
      {input.forcepos} \
      22 \
      {output.p}
    """

