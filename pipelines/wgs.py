
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

INDEXCOVA      =     expand("qc/indexcov/{group}_indexcovALL/{group}_indexcovALL-indexcov.bed.gz", group=GROUP)
BAMSAMPLE       =     expand("sample_regions/sample_{sample}.bam", sample=SAMPLES)
BAMSAMPLE2       =     expand("sample_regions/sample_chrY_{sample}.bam", sample=SAMPLES)
MOSDEPTH       =     expand("mosdepth/{sample}.mosdepth.summary.txt", sample=SAMPLES)
BIGWIGS       =     expand("mosdepth/{sample}.bigWig", sample=SAMPLES)


rule all:
    input:
        INDEXCOVA,
        BAMSAMPLE,
        BAMSAMPLE2,
        BIGWIGS

def get_sex(wildcards):
  """
  Get FNA file for bwa mapping depending on if mouse or hum an
  """
  if "hg38" in GENOME:
    return("chrX,chrY")
  else:
    return("X,Y")

rule indexcovAll:
  ###### use bed index to find genome wide coverage of all samples together
  ###### The PCA plot generated in index.html can be useful for QC
  ###### testing correct sex, and that input has high ChrM coverage
  input:
    b=expand(join(FASTQ_DIR, "{sample}.bam"), sample=SAMPLES),
    idx=expand(join(FASTQ_DIR, "{sample}.bam.bai"), sample=SAMPLES),
  output:
    "qc/indexcov/{GROUP}_indexcovALL/{GROUP}_indexcovALL-indexcov.bed.gz"
  params:
    d="qc/indexcov/{GROUP}_indexcovALL",
    sex=get_sex
  shell:
    "goleft indexcov --sex {params.sex} -d {params.d} {input.b}"

rule sampleBed:
  output:
    "wgs_sample.bed",
    "chrY.bed"
  threads:
    1
  shell:
    """
echo -e 'chr5\t137631090\t137979329
chr9\t117255004\t118634993
chr17\t15451976\t15654425
chr15\t79517885\t79648720
chr2\t119088323\t119161774' > wgs_sample.bed

echo -e 'chrY\t1\t91744698' > chrY.bed
    """

rule sampleBAM2:
  input:
    b=join(FASTQ_DIR, "{sample}.bam"),
    idx=join(FASTQ_DIR, "{sample}.bam.bai"),
    reg="chrY.bed"
  output:
    "sample_regions/sample_chrY_{sample}.bam"
  threads:
    1
  shell:
    """
    samtools view {input.b} -L {input.reg} -b -o {output}
    samtools index {output}
    """

rule sampleBAM:
  input:
    b=join(FASTQ_DIR, "{sample}.bam"),
    idx=join(FASTQ_DIR, "{sample}.bam.bai"),
    reg="wgs_sample.bed"
  output:
    "sample_regions/sample_{sample}.bam"
  threads:
    1
  shell:
    """
    samtools view {input.b} -L {input.reg} -b -o {output}
    samtools index {output}
    """


rule mosDepth:
  input:
    b=join(FASTQ_DIR, "{sample}.bam"),
    idx=join(FASTQ_DIR, "{sample}.bam.bai"),
  output:
    "mosdepth/{sample}.mosdepth.summary.txt"
  params:
    prefix="mosdepth/{sample}",
    window=500
  shell:
    """
    mosdepth -n --fast-mode --by {params.window} {params.prefix} {input.b}
    """

rule sort:
  input:
    "mosdepth/{sample}.regions.bed.gz",
  output:
    "mosdepth/{sample}.regions.sorted.bed",
  shell:
      """
      zcat {input} | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
      """

rule graphtobigwig:
  #
  # Compress bedgraph to Bigwig for IVG visualisation
  #
  input:
    graph="mosdepth/{sample}.regions.sorted.bed",
    sizes="mm10.chrom.txt" #cut -f 1,2 mm10.fa.fai > mm10.chrom.txt
  output:
    "mosdepth/{sample}.bigWig"
  threads:
    1
  shell:
    """
    bedGraphToBigWig {input.graph} {input.sizes} {output}
    """

