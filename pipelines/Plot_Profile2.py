

# To test if enrichment is due to specific lengths of fragments

#snakemake --cores 15 --snakefile pipelines/Plot_Profile2.py -npr

#snakemake --snakefile pipelines/Plot_Profile2.py --dag | dot -Tsvg > pipelines/Plot_Profile2_dag.svg

# requires:
# ../single-cell/sequencing/metadata/hg38_sizes.chrom
# FragPos/Fragment_Position_{sample}.sorted.bed

configfile: 'pipelines/config.yml'
include: "sample_pairings.py"

import itertools
from os.path import join
METADATA_DIR = config["metadata_dir"]
GENOME = 'hg38'

# NB remember to add to regionNames: in config.yml if you add to this list
locations = ["SingleBasePeaks.NA15-SRR5627138_AND_NA15-SRR5627139_vs_NA15-SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded",
              "SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded",
              "SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL",
              "SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL",
              "SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL",
              "SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108.p0.000001.sep250.ALL"]



rule all:
  input:
    #[f"bwplots/{sample}_VS_{control}_AT_{locations}.pdf" for sample, control, locations in [[samples[0],samples[1],location] for samples,location in list(itertools.product(sample_pairings, locations))]],
    #[f"bwprofilesNorm/{sample}_VS_{control}_AT_{locations}.profile" for sample, control, locations in [[samples[0],samples[1],location] for samples,location in list(itertools.product(sample_pairings, locations))]],
    [f"results/bwplots/AllSamples_AT_{location}.pdf" for location in locations]




rule removeUTH3K4:
  """
  Remove peak if it's center overlaps with region (CI peak) of H3K4me3 peak
  Filter cov_input >5
  remove top 5 by likelihood (only really nececcary for chimp due to lack of propper input contol)
  Sort by enrichment descending
  Remove peaks with greater than 99.9% input  
  """
  input:
    bed = "data/peaks/{sample}.bed",
    UTh3k4 = "data/peaks/SingleBasePeaks.NA15-SRR5627150_AND_NA15-SRR5627151_vs_NA15-SRR5627142.p0.000001.sep250.ALL.bed",
  output:
    bed = "data/QCpeaks/{sample}_QCfiltered.bed",
  shell:
    """
    # calculate 99.9% percentile of input
    THRESHOLD=$(sort -k8 -n {input} |
    awk '{{all[NR] = $0}} END{{print all[int(NR*0.999 - 0.5)]}}' |
    cut -f8)

    cut -f 1,4,5 {input.UTh3k4} | tail -n +2 |
    bedtools slop -i stdin -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 2500 |
    bedtools intersect -v -a  {input.bed} \
                          -b stdin |
    awk ' $8 >= 5 ' |
    awk -v var="$THRESHOLD" ' $8 < var ' |
    sort -k10 -n -r |
    tail -n +6 |
    sort -k9 -n -r > {output.bed}
    """


rule downloadFASTA:
  """
  To enable extraction of FASTA to infer motifs within
  """
  output:
    fa = "data/motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    fai = "data/motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai"
  shell:
    """
    # Download Fasta
    wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -O data/motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    gunzip data/motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    sed -i 's/>/>chr/g' {output.fa}
    samtools faidx {output.fa}
    """


rule extractFASTA:
  """
  For finding Prdm9 motifs in
  """
  input:
    bed = "data/QCpeaks/{sample}_QCfiltered.bed",
    fasta = "data/motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
  output:
    fasta = "data/motifs/{sample}.fasta",
  shell:
    """
    # | cut -f 1-3 # for names based on chr loc
    bedtools slop -i {input.bed} -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 150 |
    bedtools getfasta -s -fi {input.fasta} -bed stdin -name > {output.fasta}
    """


rule downloadMotifs:
  """
  Get Altemose Motifs to infer motif positions in new peak regions
  """
  output:
    "data/motifs/Human_Motif_Results_Final_iter240.r"
  shell:
    """
    wget https://github.com/altemose/PRDM9-map/raw/master/3_MotifFinding/Human_Motif_Results_Final_iter240.r -O {output}
    """


rule findMotifs:
  """
  Locate Prdm9 motif within peak, and recenter peak there (& strand according to motif orientation)
  """
  input:
    fasta = "data/motifs/{sample}.fasta",
    motifs = "data/motifs/Human_Motif_Results_Final_iter240.r",
    bed = "data/QCpeaks/{sample}_QCfiltered.bed",
  output:
    bed = "data/QCpeaks/{sample}_MotifCenteredStranded.bed",
    plots = "data/motifs/P9_motif_locations_{sample}.pdf"
  shell:
    """
    Rscript pipelines/center_by_motif.R {input.fasta} {input.motifs} {input.bed} {output.bed} {output.plots}
    """


rule quantiseBeds:
  """
  Split peaks into quartiles
  if not Human Prdm9 allele, skip motif finding and centering/stranding, use the QCfiltered.bed directly
  """
  input:
    lambda wc: "data/QCpeaks/{locations}.bed" if "MotifCenteredStranded" in wc.locations else "data/QCpeaks/{locations}_QCfiltered.bed" #config["prdm9_allele"][]=="Human"
  output:
    expand("data/QCpeaks/{{locations}}_Q0{quantile}.bed", quantile=[0, 1, 2, 3])
  params:
    "data/QCpeaks/{locations}_Q"
  shell:
    """
    split -d -n l/4 --additional-suffix=".bed" {input} {params}
    """


rule mergeFragPos:
  """
  Merge to fragment position bed files e.g. for replicates
  """
  input:
    a="data/FragPos/Fragment_Position_{sampleA}.sorted.bed",
    b="data/FragPos/Fragment_Position_{sampleB}.sorted.bed"
  output:
    sample = "data/FragPos/Fragment_Position_{sampleA}_AND_{sampleB}.sorted.bed",
  threads:
    5
  shell:
    """
    cat {input.a} {input.b} | LC_ALL=C sort -k1,1V -k2,2n -S5G --parallel=5 > {output}
    """

rule totalBases:
  """
  Calculate total number of bases the fragments cover - for normalisation purposes
  """
  input:
    "data/FragPos/Fragment_Position_{sample}.sorted.bed"
  output:
    "data/FragPos/Fragment_Position_{sample}.total"
  shell:
    "grep -P 'chr[0-9X]+\t' {input} | awk 'BEGIN{{size=0;}}{{size = size + $3 - $2;}}END{{print size;}}' > {output}"


rule bedtograph:
  """
  To create bedgraph after merging samples (single sample bedgraph created in Map_Reads.py snakefile)
  """
  input:
    bed="data/FragPos/Fragment_Position_{sample}.sorted.bed",
    sizes=join(METADATA_DIR, "{GENOME}_sizes.chrom".format(GENOME=GENOME))
  output:
    "data/bedgraphs/depth_{sample}.bedgraph"
  threads:
    5
  shell:
    "bedtools genomecov -bg -i {input.bed} -g {input.sizes} | LC_ALL=C sort -k1,1 -k2,2n -S5G --parallel=5 > {output}"


rule graphtobigwig:
  """
  To create bigwig from bedgraph after merging samples (single sample bedgraph created in Map_Reads.py snakefile)
  """
  input:
    graph="data/bedgraphs/depth_{sample}.bedgraph",
    sizes=join(METADATA_DIR, "{GENOME}_sizes.chrom".format(GENOME=GENOME))
  output:
    "data/bedgraphs/depth_{sample}.bigWig"
  threads:
    1
  shell:
    "bedGraphToBigWig {input.graph} {input.sizes} {output}"


rule makeBED6:
  """
  Reformat BED regions because bwtool is fussy and wont' accept anything but integers in bed position 4 and 5
  Also remove regions if there's another within Xkbp so they're not double counted / pollute the surrounding regions
  """
  input:
    "data/QCpeaks/{locations}"
  output:
    "data/bed6/{locations}"
  params:
    4000 # upstream+downstream
  shell:
    """
    sort -k1,1 -k2,2n {input} | \
    mergeBed -i stdin -d {params} -c 1,1,6 -o count,count,distinct | \
    awk '$4 < 2 {{print $0}}' > {output}
    
    #awk -v OFS='\t' '{{print $1,$2,$3,"0","0",$6; }}' {input} > {output}
    """


rule randomBED:
  """
  randomly select non overlapping regions, to check normalisation works
  """
  input:
    "../single-cell/sequencing/metadata/hg38_sizes.chrom"
  output:
    "data/other/random.bed"
  shell:
    """
    head -23 {input} |
    bedtools random -g stdin -n 10000 -l 1 -seed 72346 |
    sort -k1,1 -k2,2n | \
    mergeBed -i stdin -d 4000 -c 1,1,6 -o count,count,distinct | \
    awk '$4 < 2 {{print $0}}' > {output}
    """

# sort -k1,1 -k2,2n random.bed | bedtools merge -i stdin > random_deoverlap.bed
# does nothing


rule bigwigProfile:
  """
  calculate mean coverage ("profile") over regions
  """
  input:
    locations=lambda wc: "data/bed6/{locations}.bed" if(wc.locations=="top_transcripts_ens") else expand("data/bed6/{{locations}}_Q0{quantile}.bed", quantile=[0, 1, 2, 3]),
    random = "data/other/random.bed",
    sample = "data/bedgraphs/depth_{sample}.bigWig",
  output:
    "data/bwprofiles/{sample}_AT_{locations}.profile"
  params:
    width=2000
  threads:
    1
  run:
    #"bwtool aggregate {params.width}:{params.width} {params.bedlist} {input.sample} {output} -fill=0 -firstbase"
    
    # because space sperated file list isn't accepted
    bedlist = ",".join(input.locations + [input.random])
    shell("bwtool aggregate {params.width}:{params.width} {bedlist} {input.sample} {output} -fill=0 -firstbase")


rule NormaliseProfile:
  """
  Normalise average fragment depth by an Input
  """
  input:
    sample =  "data/bwprofiles/{sample}_AT_{locations}.profile",
    sample_t =  "data/FragPos/Fragment_Position_{sample}.total",
    control = "data/bwprofiles/{control}_AT_{locations}.profile",
    control_t = "data/FragPos/Fragment_Position_{control}.total"
  output:
    "data/bwprofilesNorm/{sample}_VS_{control}_AT_{locations}.profile"
  shell:
    """
    # calculate normalisation ratio
    ts=$(cat {input.sample_t})
    tc=$(cat {input.control_t})
    ratio=$(echo "scale=3 ; $tc/$ts" | bc)
    
    # normalise sample by control + lib size, two sets of 6 columns, pos, Q1-4 & random.
    paste {input.sample} {input.control} | awk -v r="$ratio" -v OFS='\t' '{{print $1, r*$2/$8, r*$3/$9, r*$4/$10, r*$5/$11, r*$6/$12;}}' > {output}
    """


rule PlotProfile:
  """
  Plot profile graph for a single sample-input pair
  """
  input:
    sample =  "data/bwprofiles/{sample}_AT_{locations}.profile",
    sample_t =  "data/FragPos/Fragment_Position_{sample}.total",
    control = "data/bwprofiles/{control}_AT_{locations}.profile",
    control_t = "data/FragPos/Fragment_Position_{control}.total"
  output:
    "results/bwplots/{sample}_VS_{control}_AT_{locations}.pdf"
  threads:
    1
  shell:
    """
    ts=$(cat {input.sample_t})
    tc=$(cat {input.control_t})
    Rscript pipelines/plot_profile.R {output} {input.sample} $ts {input.control} $tc
    """


rule MultiProfilePlot:
  """
  Plot profile graph for all normalised samples for a given set of regions
  """
  input:
    samplePair =  [f"data/bwprofilesNorm/{sample}_VS_{control}_AT_{{locations}}.profile" for sample, control, name in sample_pairings]
  output:
    "results/bwplots/AllSamples_AT_{locations}.pdf"
  params:
    sampleName = " ".join([f"'{name}'" for s,c,name in sample_pairings]),
    regionsName =  lambda wc: config["regionNames"][wc.locations],
    width = 20,
    height = 20,
  shell:
    """
    Rscript pipelines/MultiProfilePlot.R {output} '{params.regionsName}' {params.width} {params.height} {input.samplePair} {params.sampleName}
    """

