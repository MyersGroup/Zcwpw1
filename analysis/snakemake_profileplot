
# for each of the bed files in deeptools/beds/ specified as samples in config.yml
# make a profile plot for each of the bigwig files listed below

# snakemake --cores 15 --snakefile snakemake_profileplot -npr

configfile: 'config.yml'

rule all:
  input:
    [f"deeptools/plots/{location}.pdf" for location in config["samples"]]


# remove peak if it's center overlaps with region (CI peak) of H3K4me3 peak
# cov_input >5
# sort by enrichment descending
# remove peaks with greater than 99.9% input
rule removeUTH3K4:
  input:
    bed = "peaks/{sample}.bed",
    UTh3k4 = "peaks/SingleBasePeaks.SRA_Altemose2015_SRR5627150_AND_SRA_Altemose2015_SRR5627151_vs_SRA_Altemose2015_SRR5627142.p0.000001.sep250.ALL.bed",
  output:
    bed = "motifs/{sample}_QCfiltered.bed",
  shell:
    """
    # calculate 99.9% percentile of input
    THRESHOLD=$(sort -k8 -n {input} |
    awk '{{all[NR] = $0}} END{{print all[int(NR*0.999 - 0.5)]}}' |
    cut -f8)

    cut -f 1,4,5 {input.UTh3k4} | tail -n +2 |
    bedtools slop -i stdin -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 1000 |
    bedtools intersect -v -a  {input.bed} \
                          -b stdin |
    awk ' $8 >= 5 ' |
    awk -v var="$THRESHOLD" ' $8 < var ' |
    sort -k9 -n -r > {output.bed}
    """


rule downloadFASTA:
  output:
    fa = "motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    fai = "motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fai"
  shell:
    """
    # Download Fasta
    wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -O motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    gunzip motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    sed -i 's/>/>chr/g' {output.fa}
    samtools faidx {output.fa}
    """

rule downloadMotifs:
  output:
    "motifs/Human_Motif_Results_Final_iter240.r"
  shell:
    """
    wget https://github.com/altemose/PRDM9-map/raw/master/3_MotifFinding/Human_Motif_Results_Final_iter240.r -O {output}
    """

rule extractFASTA:
  input:
    bed = "motifs/{sample}_QCfiltered.bed",
    fasta = "motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
  output:
    fasta = "motifs/{sample}.fasta",
  shell:
    """
    # | cut -f 1-3 # for names based on chr loc
    bedtools slop -i {input.bed} -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 150 |
    bedtools getfasta -s -fi {input.fasta} -bed stdin -name > {output.fasta}
    """

rule findMotifs:
  input:
    fasta = "motifs/{sample}.fasta",
    motifs = "motifs/Human_Motif_Results_Final_iter240.r",
    bed = "motifs/{sample}_QCfiltered.bed",
  output:
    bed = "deeptools/beds/{sample}_MotifCenteredStranded.bed",
    plots = "motifs/P9_motif_locations_{sample}.pdf"
  shell:
    """
    Rscript center_by_motif.R {input.fasta} {input.motifs} {input.bed} {output.bed} {output.plots}
    """

# if not Human Prdm9 allele, skip motif finding and centering/stranding, use the QCfiltered.bed directly
rule quantiseBeds:
  input:
    lambda wc: "deeptools/beds/{locations}.bed" if config["prdm9_allele"][wc.locations]=="Human" else "motifs/{locations}_QCfiltered.bed"
  output:
    expand("deeptools/beds/{{locations}}.bed0{quantile}", quantile=[0, 1, 2, 3])
  params:
    "deeptools/beds/{locations}.bed"
  shell:
    """
    split -d -n l/4 {input} {params}
    """


rule gettopTSS:
  input:
  output:
    bed="deeptools/beds/top_transcripts_ens2.bed",
    topgenes="deeptools/beds/top_HEK293_genes.txt"
  shell:
    """
    wget https://www.proteinatlas.org/download/transcript_rna_celline.tsv.zip
    unzip transcript_rna_celline.tsv.zip

    # cut HEK239 expression
    # get top 20k expressed genes
    # for each gene get top expressed transcript
    # cut ensembl transcript ID
    cut -f1-2,44 transcript_rna_celline.tsv | sort -k3rn | head -20000 | sort -u -k1,1 -k3n,1 | cut -f2 > {output.topgenes}

    wget -qO- http://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz \
      | gunzip --stdout - \
      | awk '$3 == "transcript"' - \
      | grep "protein_coding" - \
      | convert2bed -i gtf - \
      | grep -Ff {output.topgenes} - \
      > {output.bed}
    """

  #wget -qO- http://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz \
  #    | gunzip --stdout - \
  #    | awk '$3 == "gene"' - \
  #    | grep "protein_coding" - \
  #    | convert2bed -i gtf - \
  #    > deeptools/beds/genes_prot_coding_ens.bed

  # alternatives
  #ftp://ftp.ebi.ac.uk/pub/databases/gencode/ensembl_ftp_files/ens_94_human/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
  #ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
  #ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz
  # abandoned attempt at getting genes expressed in HEK293T
  # wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
  # wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99407/suppl/GSE99407%5FRNAseq%5FCuffDiff%2Egenes%2Ediff%2Egz
  #
  # cuffdiff <- fread("zcat < ~/Downloads/GSE99407_RNAseq_CuffDiff.genes.diff.gz")
  # head(unique(cuffdiff[sample_1=="Human"][order(-value_1)]$gene),250)


rule deeptoolsMatrix:
  input:
    lambda wc: "deeptools/beds/{locations}.bed" if(wc.locations=="top_transcripts_ens") else expand("deeptools/beds/{{locations}}.bed0{quantile}", quantile=[0, 1, 2, 3])
  output:
    "deeptools/matrices/{locations}.gz"
  threads:
    15
  shell:
    """
    computeMatrix reference-point \
    -S deeptools/bigwigs/WTCHG_538916_221156_vs_WTCHG_538916_217108_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_538916_223180_vs_WTCHG_538916_220144_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_538916_223180_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_538916_223180_vs_WTCHG_538916_221156_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_538916_224192_vs_WTCHG_538916_221156_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627146_AND_SRA_Altemose2015_SRR5627147_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627145_AND_SRA_Altemose2015_SRR5627144_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2_bs5.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627149_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627149_vs_SRA_Altemose2015_SRR5627143_CPM_log2_NEW.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627150_vs_SRA_Altemose2015_SRR5627142_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627148_vs_SRA_Altemose2015_SRR5627142_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_659233_256196_vs_WTCHG_659233_255184_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_659233_234122_vs_WTCHG_659233_233110_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_659233_235134_vs_WTCHG_659233_233110_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_659233_250124_vs_WTCHG_659233_249112_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_659233_251136_vs_WTCHG_659233_249112_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_659233_253160_vs_WTCHG_659233_252148_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_659233_254172_vs_WTCHG_659233_252148_CPM_log2.bw \
    -R {input} \
    -a 2000 \
    -b 2000 \
    -p 15 \
    --skipZeros \
    -o {output}
    """

rule deeptoolsProfile:
  input:
    matrix="deeptools/matrices/{locations}.gz",
  output:
    "deeptools/plots/{locations}.pdf"
  params:
    regions = lambda wc: "Q1 Q2 Q3 Q4" if(wc.locations!="top_transcripts_ens") else "Expressed_Promoters",
    ymax = lambda wc: config["ymax"][wc.locations],
    center_label = lambda wc: config["center_label"][wc.locations]
  threads:
    1
  shell:
    """
  plotProfile \
    -m {input.matrix} \
    -out {output} \
    --refPointLabel "{params.center_label}" \
    --samplesLabel ChpHA_ZHA-vs-In_ZHA_CPM ChpHA_ZHA_hP9V5-vs-Input ChpHA_ZHA_hP9V5-vs-ComboInput ChpHA_ZHA_hP9V5-vs-ChpHA_ZHA ChpHA_ZHA_cP9V5-vs-ChpHA_ZHA  hPrdm9_Nterm-vs-Input hPrdm9_Cterm-vs-Input cPrdm9_Cterm-vs-Input H3K4wPrdm9_v_in H3K4wPrdm9_v_in_bs5 H3K36wPrdm9_v_in H3K36wPrdm9_v_in_SB H3K4_untransfected_v_in H3K36_untransfected_v_in ChpV5_hP9-ZFonlyV5 ChpHA_ZHA+hP9-ZonlyV5 ChpH3K4me3_ZHA+hP9-ZonlyV5 ChpHA_ZHA+hP9V5 ChpH3K4me3_ZHA+hP9V5 ChpHA_ZHA+cP9V5 ChpH3K4me3_ZHA+cP9V5 \
    --numPlotsPerRow 4 \
    --regionsLabel {params.regions} \
    --yMax {params.ymax}
    """
