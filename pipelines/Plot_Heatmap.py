

configfile: 'pipelines/config.yml'
include: "sample_pairings.py"

import itertools
from os.path import join
METADATA_DIR = config["metadata_dir"]
GENOME = 'hg38'

locations = ["SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded",
              "SingleBasePeaks.NA15-SRR5627138_AND_NA15-SRR5627139_vs_NA15-SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded",
              "SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL",
              "SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL",
              "SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL",
              "SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108.p0.000001.sep250.ALL"]

rule all:
  input:
    [f"bwmatrices/{sample}_AT_{locations}_Q00.bwm" for sample, control, locations in [[samples[0],samples[1],location] for samples,location in list(itertools.product(sample_pairings, locations))]],
    [f"bwmatrices/{control}_AT_{locations}_Q00.bwm" for sample, control, locations in [[samples[0],samples[1],location] for samples,location in list(itertools.product(sample_pairings, locations))]]


rule bigwigMatrix:
  """
  Extract regions
  """
  input:
    locations=lambda wc: "bed6/{locations}_Q00.bed" if(wc.locations=="top_transcripts_ens") else "bed6/{locations}_Q00.bed",
    sample = "bedgraphs/depth_{sample}.bigWig",
  output:
    "bwmatrices/{sample}_AT_{locations}_Q00.bwm"
  params:
    width=2000
  threads:
    1
  shell:
    """
    bwtool matrix -fill=0 -decimals=1 -tiled-averages=5 {params.width}:{params.width} {input.locations} {input.sample} {output}
    """


rule PlotProfile:
  """
  Plot profile graph for a single sample-input pair
  """
  input:
    sample =  "bwprofiles/{sample}_AT_{locations}.profile",
    sample_t =  "FragPos/Fragment_Position_{sample}.total",
    control = "bwprofiles/{control}_AT_{locations}.profile",
    control_t = "FragPos/Fragment_Position_{control}.total"
  output:
    "bwplots/{sample}_VS_{control}_AT_{locations}.pdf"
  threads:
    1
  shell:
    """
    ts=$(cat {input.sample_t})
    tc=$(cat {input.control_t})
    Rscript pipelines/plot_profile.R {output} {input.sample} $ts {input.control} $tc
    """

