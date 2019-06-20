# snakemake --cores 15 --snakefile Snakefile_peaks -npr


# Combine Inputs, needs to be done manually atm

# cat MAPeakCaller/Fragment_Position_WTCHG_538916_217108.sorted.bed \
#     MAPeakCaller/Fragment_Position_WTCHG_538916_220144.sorted.bed | \
#     LC_ALL=C sort -k1,1V -k2,2n -S5G --parallel=5 \
#     > MAPeakCaller/Fragment_Position_WTCHG_538916_217108_AND_WTCHG_538916_220144.sorted.bed

# LC_ALL=C sort -k1,1V -k2,2n -S5G --parallel=5 \
# > Fragment_Position_${sample}.sorted.bed



from os.path import join
METADATA_DIR = '../single-cell/sequencing/metadata'
GENOME = 'hg38'

# chip, control
# when there are two biological replicates, seperate them using the special string '_AND_'
# when there are pseudoreplicates, just enter the common sample name, assumes the PR1 PR2 file ending convention was used

comparisons = [("WTCHG_538916_221156", "WTCHG_538916_217108"),  #Zcw only Chip vs Input
              ("WTCHG_538916_223180", "WTCHG_538916_220144"), #Cotransfected Chip vs Cotransfected Input
              ("WTCHG_538916_223180", "WTCHG_538916_217108_AND_WTCHG_538916_220144"), #Cotransfected Chip vs Combined Input

              ("WTCHG_538916_223180", "WTCHG_538916_221156"),  #Chip with Prdm9 vs Chip
              ("WTCHG_538916_224192", "WTCHG_538916_221156"), # chimpPrdm9+Zcw vs chip
              ("WTCHG_538916_224192", "WTCHG_538916_217108_AND_WTCHG_538916_220144"), # chimpPrdm9+Zcw vs Combined hInput

              # Human
              ("SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139", "SRA_Altemose2015_SRR5627140"), # hPrdm9-Nterm vs YFP_HumanPRDM9.Input

              ("SRA_Altemose2015_SRR5627141", "SRA_Altemose2015_SRR5627140"), # YFP_HumanPRDM9.antiH3K4me3 vs YFP_HumanPRDM9.Input
              #("SRA_Altemose2015_SRR5627141", "SRA_Altemose2015_SRR5627136_AND_SRA_Altemose2015_SRR5627135"), # YFP_HumanPRDM9.antiH3K4me3 vs Untransfected.antiH3K4me3

              ("SRA_Altemose2015_SRR5627146", "SRA_Altemose2015_SRR5627143"), # hPrdm9 HA vs Input
              ("SRA_Altemose2015_SRR5627147", "SRA_Altemose2015_SRR5627143"), # hPrdm9 V5 vs Input (HA)
              ("SRA_Altemose2015_SRR5627146_AND_SRA_Altemose2015_SRR5627147", "SRA_Altemose2015_SRR5627143"), # hPrdm9-Cterm HA & V5 vs Input (HA)

              # Chimp
              ("SRA_Altemose2015_SRR5627145", "SRA_Altemose2015_SRR5627143"), # cPrdm9 V5 vs Input (HA)
              ("SRA_Altemose2015_SRR5627144", "SRA_Altemose2015_SRR5627143"), # cPrdm9 HA vs (h) Input
              ("SRA_Altemose2015_SRR5627145_AND_SRA_Altemose2015_SRR5627144", "SRA_Altemose2015_SRR5627143"), # cPrdm9 V5 & HA vs Input (HA)

              # Histone Mods
              ("SRA_Altemose2015_SRR5627152_AND_SRA_Altemose2015_SRR5627153", "SRA_Altemose2015_SRR5627143"), # H3K4 peaks wPrdm9
              ("SRA_Altemose2015_SRR5627149", "SRA_Altemose2015_SRR5627143"), # H3K36 peaks wPprdm9
              ("SRA_Altemose2015_SRR5627150_AND_SRA_Altemose2015_SRR5627151", "SRA_Altemose2015_SRR5627142"), # untransfected H3K4 peaks
              ("SRA_Altemose2015_SRR5627148", "SRA_Altemose2015_SRR5627142"), # untransfected H3K36 peaks


              # Chip H3K4 & ZF only
              ("WTCHG_659233_256196", "WTCHG_659233_255184"), # ChIPV5_hPRDM9-ZFonly-V5 vs input

              ("WTCHG_659233_253160", "WTCHG_659233_252148"), # ChIPHA_ZCWPW1-HA+cPRDM9-V5 vs input
              ("WTCHG_659233_254172", "WTCHG_659233_252148"), # ChIPH3K4me3_ZCWPW1-HA+cPRDM9-V5 vs input

              ("WTCHG_659233_250124", "WTCHG_659233_249112"), # ChIPHA_ZCWPW1-HA+hPRDM9-V5 vs input
              ("WTCHG_659233_251136", "WTCHG_659233_249112"), # ChIPH3K4me3_ZCWPW1-HA+hPRDM9-V5 vs input

              ("WTCHG_659233_234122", "WTCHG_659233_233110"), # ChIPHA_ZCWPW1-HA+hPRDM9-ZFonly-V5 vs input
              ("WTCHG_659233_235134", "WTCHG_659233_233110")] # ChIPH3K4me3_ZCWPW1-HA+hPRDM9-ZFonly-V5 vs input

rule all:
  input:
    ["peaks/SingleBasePeaks.{chip}_vs_{control}.p0.000001.sep250.ALL.bed".format(
        chip = chip_id,
        control=control_id) for chip_id, control_id in comparisons],
    ["deeptools/bigwigs/{chip}_vs_{control}_CPM_log2.bw".format(
        chip = chip_id,
        control=control_id) for chip_id, control_id in comparisons]
    #["deeptools/bigwigs/SBPsample/{chip}_vs_{control}_CPM_log2.bw".format(
    #    chip = chip_id,
    #    control=control_id) for chip_id, control_id in comparisons]


def get_chips(wildcards):
  """
  When the input is required to be two file, return either the two biological replicates, or the two pseudoreplicates
  """
  if "_AND_" in wildcards.chip:
    two_samples = wildcards.chip.split("_AND_")
    return(list([f'MAPeakCaller/Fragment_Position_{two_samples[0]}.sorted.bed',
                f'MAPeakCaller/Fragment_Position_{two_samples[1]}.sorted.bed']))
  else:
    return(list([f'MAPeakCaller/Fragment_Position_{wildcards.chip}.sorted.bed.PR1',
                f'MAPeakCaller/Fragment_Position_{wildcards.chip}.sorted.bed.PR2']))

def split_chip_bam(wildcards):
  """
  When the input is required to be two file, split the wildcard on AND
  """
  if "_AND_" in wildcards.chip:
    two_samples = wildcards.chip.split("_AND_")
    return(list([f'filtered/{two_samples[0]}.bam',
                f'filtered/{two_samples[1]}.bam']))
  else:
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


#rule call_peaks:
#  input:
#    chips = get_chips ,
#    control="MAPeakCaller/Fragment_Position_{control}.sorted.bed",
#    chrsizes=join(METADATA_DIR, "hg38_sizes.chrom")
#  output:
#    c="peaks/Constants.{chip}_vs_{control}.tsv",
#    p="peaks/SingleBasePeaks.{chip}_vs_{control}.p0.000001.sep250.ALL.bed"
#  threads:
#    15
#  shell:
#    """
#    sh MAPeakCaller.sh \
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
    control="MAPeakCaller/Fragment_Position_{control}.sorted.bed",
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
    chips = split_chip_bam ,
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


