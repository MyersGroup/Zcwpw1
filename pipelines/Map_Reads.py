# snakemake --cores 15 -npr --config GROUP="659233"
# snakemake --cores 15 -npr --config GROUP="538916"
# snakemake --cores 15 -npr --config GROUP="594404"
# snakemake --cores 15 -npr --config GROUP="Altemose2015"
# snakemake --cores 15 -npr --config GROUP="Dmc1_r1"

# snakemake --dag | dot | display

from os.path import join

configfile: 'pipelines/config.yml'

METADATA_DIR = config["metadata_dir"]

GROUP = config["GROUP"]

FASTQ_DIR = config["fastq_dir"][GROUP]
print(FASTQ_DIR)
GENOME = config["genome"][GROUP]

SAMPLES, PAIRS = glob_wildcards(join(FASTQ_DIR, "{samples}_{pairs}.fastq.gz"))

SAMPLES = list(set(SAMPLES))
print(SAMPLES)

##################
# Desired output
##################

FASTQC_REPORTS  =     expand("qc/fastqc/{sample}_{pair}_fastqc.html", pair=[1,2], sample=SAMPLES)
LIBCOMPLEX      =     expand("qc/LibComplex/{sample}_LibComplex.txt", sample=SAMPLES)
BIGWIG          =     expand("bedgraphs/depth_{sample}.bigWig", sample=SAMPLES)
FRAGPOS         =     expand("MAPeakCaller/Fragment_Position_{sample}.sorted.bed.PR{replicate}", replicate=[1,2], sample=SAMPLES)
FLAGSTAT        =     expand("qc/FlagStat/{sample}_FlagStat.txt", sample=SAMPLES)
FLAGSTAT2       =     expand("qc/FlagStat/{sample}_preQC_FlagStat.txt", sample=SAMPLES)
BLACKLISTS      =     expand("blacklists/{sample}_blacklist.bed", sample=SAMPLES)
MOSDEPTH1       =     expand("qc/mosdepth/{sample}.mosdepth.global.dist.txt", sample=SAMPLES)
MOSDEPTH2       =     expand("qc/mosdepth/raw_{sample}.mosdepth.global.dist.txt", sample=SAMPLES)
MOSDEPTHP       =     expand("qc/mosdepth/mosdepth_{group}.pdf", group=GROUP)
FRAGSIZE        =     expand("qc/deeptools/{group}_frag_size.pdf", group=GROUP)
FINGERPRINT     =     expand("qc/deeptools/{group}_fingerprint.pdf", group=GROUP)
INDEXCOVA      =     expand("qc/indexcov/{group}_indexcovALL/{group}_indexcovALL-indexcov.bed.gz", group=GROUP)
BWCPM           =     expand("deeptools/bigwigs/{sample}_CPM.bw", sample=SAMPLES)
BAMSAMPLE       =     expand("sample_regions/sample_{sample}.bam", sample=SAMPLES)
PCA             =     expand("qc/deeptools/PCA_{group}.png", group=GROUP)
COR             =     expand("qc/deeptools/Cor_{group}.png", group=GROUP)


rule all:
    input:
        FASTQC_REPORTS,
        LIBCOMPLEX,
        BIGWIG,
        FRAGPOS,
        FLAGSTAT,
        FLAGSTAT2,
        BLACKLISTS,
        MOSDEPTH1,
        MOSDEPTH2,
        #MOSDEPTHP,
        FRAGSIZE,
        FINGERPRINT,
        INDEXCOVA,
        BWCPM,
        BAMSAMPLE,
        PCA,
        COR


    message: "ChIP-seq pipeline succesfully run."

join(FASTQ_DIR, "{sample}_1.fastq.gz"),

rule fastqc:
  input:
    fq1=join(FASTQ_DIR, "{sample}_1.fastq.gz"),
    fq2=join(FASTQ_DIR, "{sample}_2.fastq.gz")
  output:
    qc1="qc/fastqc/{sample}_1_fastqc.html",
    qc2="qc/fastqc/{sample}_2_fastqc.html"
  shell:
    """
    fastqc {input.fq1} -o qc/fastqc
    fastqc {input.fq2} -o qc/fastqc
    """

# https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

rule download_hg38:
  output:
    fna=join(METADATA_DIR, "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"),
    ann=join(METADATA_DIR, "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz.ann")
  shell:
    """
    wget --directory-prefix {METADATA_DIR} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
    bwa index {METADATA_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
    """


# https://www.biostars.org/p/342482/
rule download_mm10:
  output:
    fna=join(METADATA_DIR, "Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"),
    ann=join(METADATA_DIR, "Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz.ann")
  shell:
    """
    wget --directory-prefix {METADATA_DIR} ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
    bwa index {METADATA_DIR}/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
    """

rule chr_sizes:
  input:
    join(METADATA_DIR, "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz.ann")
  output:
    join(METADATA_DIR, "hg38_sizes.chrom")
  shell:
    "awk '{{ print $2; }}' {input} | tail -n +2 | paste - - > {output}"

rule chr_sizesMm10:
  input:
    join(METADATA_DIR, "Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz.ann")
  output:
    join(METADATA_DIR, "mm10_sizes.chrom")
  shell:
    "awk '{{ print $2; }}' {input} | tail -n +2 | paste - - > {output}"


def get_fna(wildcards):
  """
  Get FNA file for bwa mapping depending on if mouse or hum an
  """
  if "hg38" in GENOME:
    return(join(METADATA_DIR, "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"))
  else:
    return(join(METADATA_DIR, "Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"))

rule bwa_map:
  input:
    g=get_fna ,
    fq=expand(join(FASTQ_DIR, "{{sample}}_{pair}.fastq.gz"), pair=[1,2])
  output:
    "bams/{sample}.bam"
  threads:
    15
  shell:
    "bwa mem -t 15 {input.g} {input.fq} | \
    samtools sort -@15 -O BAM -o {output}"

rule remove_duplicates:
  input:
    "bams/{sample}.bam"
  output:
    b="undup/{sample}.bam",
    log="undup/{sample}_metrics.txt"
  threads:
    15
  shell:
    "picard MarkDuplicates \
    INPUT={input} \
    OUTPUT={output.b} \
    METRICS_FILE={output.log} \
    REMOVE_DUPLICATES=true"


rule remove_blacklisted:
  input:
    bam="undup/{sample}.bam",
    black=join(METADATA_DIR, "{GENOME}.blacklist.bed".format(GENOME=GENOME))
  output:
    temp("tmp/{sample}_blacklisted.bam")
  shell:
    "bedtools intersect -v -a {input.bam} -b {input.black} > {output}"


rule filter_lowq:
  # https://broadinstitute.github.io/picard/explain-flags.html
  # Remove  unmapped, mate unmapped
  # not primary alignment, reads failing platform
  # Remove low MAPQ reads
  # namesorted
  input:
    "tmp/{sample}_blacklisted.bam"
  output:
    temp("tmp/{sample}_filtered_name_sorted.bam")
  threads:
    5
  shell:
    """
    samtools view -q 30 -F 3852 -f 2 -u {input} | \
      samtools sort -@ 5 -n - -o {output}
    """
#1804 + rm supp


rule fixmates:
  input:
    "tmp/{sample}_filtered_name_sorted.bam"
  output:
    temp("tmp/{sample}_filtered_name_sorted_fixmate.bam")
  threads:
    5
  shell:
    """
    samtools fixmate -@ 5 -r {input} {output}
    """

rule filter_orphans:
  # Only keep properly paired reads
  # position sorted
  input:
    "tmp/{sample}_filtered_name_sorted_fixmate.bam"
  output:
    "filtered/{sample}_namesorted.bam"
  threads:
    5
  shell:
    """
    samtools view -q 30 -F 3852 -f 2 {input} -o {output}
    """

rule position_sort:
  input:
    "filtered/{sample}_namesorted.bam"
  output:
    "filtered/{sample}.bam"
  threads:
    5
  shell:
    "samtools sort -@ 5 {input} -o {output}"


rule libcomplexity:
  input:
    "filtered/{sample}.bam"
  output:
    "qc/LibComplex/{sample}_LibComplex.txt"
  threads:
    15
  shell:
    "picard EstimateLibraryComplexity INPUT={input} OUTPUT={output}"


rule flagstat:
  input:
    "filtered/{sample}.bam"
  output:
    "qc/FlagStat/{sample}_FlagStat.txt"
  threads:
    5
  shell:
    "samtools flagstat -@ 5 {input} > {output}"

rule flagstat2:
  input:
    "bams/{sample}.bam"
  output:
    "qc/FlagStat/{sample}_preQC_FlagStat.txt"
  threads:
    5
  shell:
    "samtools flagstat -@ 5 {input} > {output}"


rule index_bam:
  input:
    "{sample}.bam"
  output:
    "{sample}.bam.bai"
  threads:
    5
  shell:
    "samtools index -@ 5 {input}"

def get_sex(wildcards):
  """
  Get FNA file for bwa mapping depending on if mouse or hum an
  """
  if "hg38" in GENOME:
    return("chrX,chrY")
  else:
    return("X,Y")


rule indexcov:
  ###### use bed index to find genome wide coverage
  ###### Used to generate blacklist regions using rule:generate_blacklist
  input:
    bam="filtered/{sample}.bam",
    idx="filtered/{sample}.bam.bai"
  output:
    "qc/indexcov/{sample}_indexcov/{sample}_indexcov-indexcov.bed.gz"
  params:
    d="qc/indexcov/{sample}_indexcov",
    sex=get_sex
  shell:
    "goleft indexcov --sex {params.sex} -d {params.d} {input.bam}"


rule indexcovAll:
  ###### use bed index to find genome wide coverage of all samples together
  ###### The PCA plot generated in index.html can be useful for QC
  ###### testing correct sex, and that input has high ChrM coverage
  input:
    b=expand("filtered/{sample}.bam", sample=SAMPLES),
    idx=expand("filtered/{sample}.bam.bai", sample=SAMPLES),
  output:
    "qc/indexcov/{GROUP}_indexcovALL/{GROUP}_indexcovALL-indexcov.bed.gz"
  params:
    d="qc/indexcov/{GROUP}_indexcovALL",
    sex=get_sex
  shell:
    "goleft indexcov --sex {params.sex} -d {params.d} {input.b}"


rule generate_blacklist:
  input:
    "qc/indexcov/{sample}_indexcov/{sample}_indexcov-indexcov.bed.gz"
  output:
    "blacklists/{sample}_blacklist.bed"
  shell:
    "Rscript gen_blacklist.R {input} {output}"


  # combinding & filtering on blacklist - no longer used on bams
  #cat blacklists/${sample}_blacklist.bed >> blacklists/${expid}_combined_blacklist.bed
  #sort -k1,1V -k2,2n -u blacklists/${expid}_combined_blacklist.bed -o blacklists/${expid}_combined_blacklist.bed
  #bedtools intersect -v -a undup/${sample}.bam -b blacklists/${expid}_combined_blacklist.bed > higcovfiltered/${sample}.bam
  #samtools index higcovfiltered/${sample}.bam


rule mosdepth:
  #
  # Compare depth/coverage distribution for raw & deduped
  #
  input:
    bam="filtered/{sample}.bam",
    idx="filtered/{sample}.bam.bai"
  output:
    "qc/mosdepth/{sample}.mosdepth.global.dist.txt"
  shell:
    "MOSDEPTH_PRECISION=5 mosdepth -n -t 4 qc/mosdepth/{wildcards.sample} {input.bam}"


rule mosdepth2:
  input:
    bam="bams/{sample}.bam",
    idx="bams/{sample}.bam.bai"
  output:
    "qc/mosdepth/raw_{sample}.mosdepth.global.dist.txt"
  shell:
    "MOSDEPTH_PRECISION=5 mosdepth -n -t 4 qc/mosdepth/raw_{wildcards.sample} {input.bam}"


rule mosdepthplot:
  input:
    expand("qc/mosdepth/raw_{sample}.mosdepth.global.dist.txt", sample=SAMPLES),
    expand("qc/mosdepth/{sample}.mosdepth.global.dist.txt", sample=SAMPLES)
  output:
    "qc/mosdepth/mosdepth_{GROUP}.pdf"
  params:
    g=GROUP
  shell:
    "Rscript plot_mosdepth.R {params}"

#Removed 479 rows containing missing values (geom_path).
#Warning message:
#In eval(jsub, SDenv, parent.frame()) : NAs introduced by coercion
#Error in mutate_impl(.data, dots) :
#  Evaluation error: `as_dictionary()` is defunct as of rlang 0.3.0.
#Please use `as_data_pronoun()` instead.
#Calls: +.gg ... as.data.frame -> mutate -> mutate.tbl_df -> mutate_impl


rule bamtobed2:
  #
  # Bam to fragment position bed file for input into peak caller
  # Deprecated as uses read length to calculate fragment/insert size (wrong),
  # & creates extra fragments (duplicates) when start of read 1 == start(coordinate)/end(real) of read 2
  #
  input:
    "filtered/{sample}.bam"
  output:
    "MAPeakCaller/Fragment_Position_{sample}.sorted.bed_DEPRECATED"
  threads:
    2
  shell:
    "samtools view {input} | perl MAPeakCaller/GetFragmentPositions.pl {output} {READ_LENGTH}"


rule bamtobed:
  #
  # Bam to fragment position bed file for input into peak caller
  #
  input:
    "filtered/{sample}_namesorted.bam"
  output:
    "MAPeakCaller/Fragment_Position_{sample}.sorted.bed"
  threads:
    5
  shell:
    """
    bedtools bamtobed -i {input} -bedpe | awk '{{print $1, $2, $6}}' OFS="\t" | LC_ALL=C sort -k1,1V -k2,2n -S5G --parallel=5 > {output}
    """

rule bedtograph:
  #
  # Create converage for IGV visualisation
  #
  input:
    bed="MAPeakCaller/Fragment_Position_{sample}.sorted.bed",
    sizes=join(METADATA_DIR, "{GENOME}_sizes.chrom".format(GENOME=GENOME))
  output:
    "bedgraphs/depth_{sample}.bedgraph"
  threads:
    5
  shell:
    "bedtools genomecov -bg -i {input.bed} -g {input.sizes} | LC_ALL=C sort -k1,1 -k2,2n -S5G --parallel=5 > {output}"


rule graphtobigwig:
  #
  # Compress bedgraph to Bigwig for IVG visualisation
  #
  input:
    graph="bedgraphs/depth_{sample}.bedgraph",
    sizes=join(METADATA_DIR, "{GENOME}_sizes.chrom".format(GENOME=GENOME))
  output:
    "bedgraphs/depth_{sample}.bigWig"
  threads:
    1
  shell:
    "bedGraphToBigWig {input.graph} {input.sizes} {output}"


rule multiBigWigSummary:
  input:
    b=expand("bedgraphs/depth_{sample}.bigWig", sample=SAMPLES),
    black=join(METADATA_DIR, "{GENOME}.blacklist.bed".format(GENOME=GENOME))
  output:
    "qc/deeptools/{GROUP}_multibigwigSummary.npz"
  threads:
    15
  shell:
    "multiBigwigSummary bins -p 15 -b {input.b} -bl {input.black} -o {output}"


rule PCA:
  # For QC
  input:
    "qc/deeptools/{GROUP}_multibigwigSummary.npz"
  output:
    p="qc/deeptools/PCA_{GROUP}.png",
    t="qc/deeptools/PCA_{GROUP}.tsv",
    p2="qc/deeptools/PCA_{GROUP}_T.png",
    t2="qc/deeptools/PCA_{GROUP}_T.tsv"
  shell:
    """
    plotPCA -in {input} -o {output.p} --outFileNameData {output.t}
    plotPCA -in {input} --transpose -o {output.p2} --outFileNameData {output.t2}
    """


rule Corr:
  # For QC
  input:
    "qc/deeptools/{GROUP}_multibigwigSummary.npz"
  output:
    p="qc/deeptools/Cor_{GROUP}.png",
    p2="qc/deeptools/Cor_{GROUP}_pear.png",
    p3="qc/deeptools/Cor_{GROUP}_pear_scatter.png"
  shell:
    """
    plotCorrelation -in {input} --corMethod spearman --skipZeros --whatToPlot heatmap -o {output.p}
    plotCorrelation -in {input} --corMethod pearson --skipZeros --whatToPlot heatmap -o {output.p2}
    plotCorrelation -in {input} --corMethod pearson --skipZeros --whatToPlot scatterplot -o {output.p3}
    """


rule pseudoreplicate:
  #
  # Create psuedoreplicate as required for peak caller if two real replicates not avaliable
  #
  input:
    "MAPeakCaller/Fragment_Position_{sample}.sorted.bed"
  output:
    expand("MAPeakCaller/Fragment_Position_{{sample}}.sorted.bed.PR{replicate}", replicate=[1,2])
  shell:
    """
    awk '{{print > (rand()<0.5 ? (FILENAME".PR1") : (FILENAME".PR2"))}}' {input}
    """


rule downloadBlacklist:
  output:
    join(METADATA_DIR,"hg38.blacklist.bed"),
    join(METADATA_DIR,"mm10.blacklist.bed")
  shell:
    """
    wget -P {METADATA_DIR} http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
    wget -P {METADATA_DIR} http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
    gzip -d {METADATA_DIR}/hg38.blacklist.bed.gz
    gzip -d {METADATA_DIR}/mm10.blacklist.bed.gz
    """

rule fingerprint:
  input:
    b=expand("filtered/{sample}.bam", sample=SAMPLES),
    idx=expand("filtered/{sample}.bam.bai", sample=SAMPLES),
    black=join(METADATA_DIR, "{GENOME}.blacklist.bed".format(GENOME=GENOME))
  output:
    "qc/deeptools/{GROUP}_fingerprint.pdf"
  threads:
    15
  shell:
    "plotFingerprint -b {input.b} --smartLabels -p 15 -plot {output} --blackListFileName {input.black}"


rule fragsize:
  input:
    b=expand("filtered/{sample}.bam", sample=SAMPLES),
    idx=expand("filtered/{sample}.bam.bai", sample=SAMPLES),
    black=join(METADATA_DIR, "{GENOME}.blacklist.bed".format(GENOME=GENOME))
  output:
    "qc/deeptools/{GROUP}_frag_size.pdf"
  threads:
    15
  shell:
    "bamPEFragmentSize -b {input.b} --maxFragmentLength 1000 --samplesLabel {input.b} -p 15 -hist {output} --blackListFileName {input.black}"


rule bamCoverageCPM:
  input:
    b="filtered/{sample}.bam",
    idx="filtered/{sample}.bam.bai"
  output:
    "deeptools/bigwigs/{sample}_CPM.bw"
  threads:
    5
  shell:
    """
    bamCoverage \
      --bam {input.b} \
      -o {output} \
      --normalizeUsing CPM \
      --centerReads \
      --extendReads \
      -p 5
    """


rule sampleBAM:
  input:
    b="filtered/{sample}.bam",
    idx="filtered/{sample}.bam.bai",
    reg="chr19_sample.bed"
  output:
    "sample_regions/sample_{sample}.bam"
  threads:
    1
  shell:
    """
    samtools view {input.b} -L {input.reg} -b -o {output}
    samtools index {output}
    """


# snakemake -npr --rulegraph | dot -Tpdf > dag.pdf
# snakemake -npr --dag --forceall | dot -Tpdf > dag.pdf
