
from os.path import join

configfile: 'pipelines/config.yml'

METADATA_DIR = config["metadata_dir"]

rule all:
    input:
        "results/PeakCallingMA.md"



rule CpG_islands:
  ####################################
  ###### Count CpG & CpG Islands per 100bp window
  ####################################
  input:
    windows = "data/forcepeaks/genome.windows.100wide.100slide.bed",
    cpgs = "data/CpG/CpG_hg38_cut.bed"
  output:
    cpg_txt = "data/CpG/cpgIslandExtUnmasked.txt",
    cpg_bed = "data/CpG/cpgIslandExtUnmasked.bed",
    cpg_bed_windows = "data/CpG/cpgIslandExtUnmasked_100bpwindows.bed",
    cpg_bed_windows_f1 = "data/CpG/cpgIslandExtUnmasked_100bpwindows_f1.bed",
    cpg_windows = "data/CpG/CpG_hg38_100bpwindows.bed"
  shell:
    """
    bedtools intersect -a {input.windows} -b {input.cpgs} -c > {output.cpg_windows}

    ##### CpG Islands per 100bp

    wget -P data/CpG/ http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz
    gzip -d data/CpG/cpgIslandExtUnmasked.txt.gz
    cut -f 2-4 {output.cpg_txt} > {output.cpg_bed}

    bedtools intersect -a {input.windows} -b {output.cpg_bed} -f 0.1 -c > {output.cpg_bed_windows}

    bedtools intersect -a {input.windows} -b {output.cpg_bed} -f 1 -c > {output.cpg_bed_windows_f1}
    """


rule human_dmc1:
  output:
    supp_file = "data/pratto_dmc1/GSE59836_Peak_data_Supplementary_File_1.txt",
    hg19 = "data/pratto_dmc1/Pratto_human_DSB_Aintersect_autosomal_hg19.bed",
    hg38 = "data/pratto_dmc1/Pratto_human_DSB_Aintersect_autosomal_hg38.bed",
    non_hg38 = "data/pratto_dmc1/Pratto_human_DSB_Aintersect_autosomal_Nonhg38.bed",
    chain = "data/genome/hg19ToHg38.over.chain.gz"
  shell:
    """
    wget -P data/genome/ https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
    wget -P data/pratto_dmc1/ ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59836/suppl/GSE59836%5FPeak%5Fdata%5FSupplementary%5FFile%5F1%2Etxt%2Egz
    gunzip {output.supp_file}.gz

    # subset to only AorB allele peaks observed in all samples & compute average strength
    awk -v OFS='\t' '{{ if ($16 == 1 && ($3-$2)<3000 && $1!="chrY" && $1!="chrX") {{ print $1, $2, $3, ($4+$5)/2, ($6+$7)/4 }} }}' \
      {output.supp_file} \
      > {output.hg19}
    # 18342
    # 18326

    liftOver {output.hg19} {output.chain} \
      {output.hg38} \
      {output.non_hg38}

    # awk -v OFS='\t' '{{ if ($11 == 1 && ($3-$2)<3000 && $1!="chrY" && $1!="chrX") {{ print $1, $2, $3, ($4+$5)/2, ($6+$7)/4 }} }}' data/dmc1/GSE59836_Peak_data_Supplementary_File_1.txt \ > data/dmc1/Pratto_human_DSB_AB_autosomal_hg19.bed
    # wc -l Pratto_human_DSB_AB_autosomal_hg19.bed
    # #20336
    #
    # awk -v OFS='\t' '{{ if ($9 == 1 && $10 == 1 && ($3-$2)<3000 && $1!="chrY" && $1!="chrX") {{ print $1, $2, $3, ($4+$5)/2, ($6+$7)/4 }} }}' data/dmc1/GSE59836_Peak_data_Supplementary_File_1.txt \ > data/dmc1/Pratto_human_DSB_AA_autosomal_hg19.bed
    # #21699
    """

rule repeats:
  input:
    hg38_fa = "data/motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
  output:
    bed = "data/repeats/Repeat_MaskerSorted.bed",
    alu = "data/repeats/Repeat_MaskerSorted_Alu.bed",
    alu_qc = "data/repeats/Repeat_MaskerSorted_Alu_QC.bed",
    not_alu = "data/repeats/Repeat_MaskerSorted_NotAlu.bed",
    alu_fa = "data/repeats/AluRepeats.fa"
  shell:
    """
    curl 'http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=738896965_zTgYMjAIHdEHcvRZcCteap7Zqbpt&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_rmsk&hgta_ctDesc=table+browser+query+on+rmsk&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED' \
      -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.14; rv:67.0) Gecko/20100101 Firefox/67.0' \
      -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
      -H 'Accept-Language: en-GB,en;q=0.5' \
      --compressed -H 'Connection: keep-alive' \
      -H 'Referer: http://genome.ucsc.edu/cgi-bin/hgTables' \
      -H 'Cookie: hguid.genome-euro=362527140_Mi0jIpneXLkgeplCtz5aa6JO13ha; hguid=673678055_VvqPPBb5iFgOJEF3WMG7M3t7iual' \
      -H 'Upgrade-Insecure-Requests: 1' > data/repeats/Repeat_Masker.bed.gz

    gunzip -d data/repeats/Repeat_Masker.bed.gz

    grep -P 'chr[0-9XY]+\t' data/repeats/Repeat_Masker.bed | sort -k1,1 -k2,2n  > {output.bed}

    # split into Alu and not
    grep 'Alu[a-zA-Z]' {output.bed} > {output.alu}
    grep -v 'Alu[a-zA-Z]' {output.bed} > {output.not_alu}

    # check Alu with no family name not present
    grep -P 'Alu\t' {output.alu} | head

    awk '$3-$2 > 250 && $3-$2 < 350' {output.alu} > {output.alu_qc}

    bedtools getfasta -s -fi {input.hg38_fa} -bed {output.alu_qc} -name > {output.alu_fa}
    """

rule per100bp:
  input:
    windows = "data/forcepeaks/genome.windows.100wide.100slide.bed",
    aluqc = "data/repeats/Repeat_MaskerSorted_Alu_QC.bed",
    notalu = "data/repeats/Repeat_MaskerSorted_NotAlu.bed",
    zcw_peaks = "data/peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.bed",
    zcw_random = "peaks/Zcw_random.bed"
  output:
    alu0p1 = "data/repeats/Repeat_MaskerSorted_Alu_QC_100bp.bed",
    alu1 = "data/repeats/Repeat_MaskerSorted_Alu_QC_100bp_f0.5.bed",
    notalu = "data/repeats/Repeat_MaskerSorted_NotAlu_100bp.bed",
    zcw_peaks = "data/peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.100bpWindows.bed",
    zcw_random_1bp = "data/peaks/Zcw_random_1bp.bed",
    zcw_random_100bp = "data/peaks/Zcw_random_1bp.100bpWindows.bed"
  shell:
    """
    # Count Alus per 100bp windows

    # f= fraction of A that must overlap

    bedtools intersect -a {input.windows} -b {input.aluqc} -f 0.1 -c > {output.alu0p1}

    bedtools intersect -a {input.windows} -b {input.aluqc} -f 1 -c > {output.alu1}

    bedtools intersect -a {input.windows} -b {input.notalu} -f 0.1 -c > {output.notalu}

    # Count Zcw peaks per 100bp (max 1 as min sep=250)
    bedtools intersect -a {input.windows} -b {input.zcw_peaks} -c > {output.zcw_peaks}

    # random peaks (made above) per 100bp
    awk -v OFS='\t' '{{print $1, $2, $2+1}}' {input.zcw_random} > {output.zcw_random_1bp}

    bedtools intersect -a {input.windows} -b {output.zcw_random_1bp} -c > {output.zcw_random_100bp}
    """

rule random_zcw:
  # Random Zcw Peaks
  input:
    chrom = "../single-cell/sequencing/metadata/hg38_sizes.chrom"
  output:
    bed = "data/peaks/Zcw_random.bed",
    sizes = "data/genome/hg38_sizes23.bed"
  shell:
    """
    head -23 {input.chrom} > {output.sizes}
    bedtools random -n 1264519 -l 300 -g {output.sizes} -seed 71346 | sort -k1,1 -k2,2n > {output.bed}
    """

rule methylation:
  ####################################
  ###### Bisulphite sequencing data
  ####################################
  input:
    windows = "data/forcepeaks/genome.windows.100wide.100slide.bed"
    chain = "data/genome/hg19ToHg38.over.chain.gz"
    cgi = "data/CpG/cpgIslandExtUnmasked.bed"
    chrom = "../single-cell/sequencing/metadata/hg38_sizes.chrom"
  output:
    cgimeth = "data/CpG/cpgIsland_Meth.bed"
    windows_namesort = "data/forcepeaks/genome.windows.100wide.100slide.namesort.bed"
    cpg_hg38 = "data/CpG/GSM1254259_HEK293-CT.CpG.hg38.bed"
    cpg_hg38_sort = "data/CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.bed"
    mean = "data/CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.bed"
    mean_chr1 = "data/CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.chr1.bed"
    bigwig = "data/CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.chr1.bigWig"
    100bp = "data/CpG/cpg_Meth_100bpwindows.bed"
  shell:
    """
    #https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000533
    #https://www.frontiersin.org/articles/10.3389/fgene.2015.00339/full
    #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1254259 (GSE51867)
    wget -P data/CpG/ ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1254nnn/GSM1254259/suppl/GSM1254259%5FHEK293%2DCT%2Ecout%2Etxt%2Egz

    # only CpG meth
    zcat data/CpG/GSM1254259_HEK293-CT.cout.txt.gz | awk ' $4 == "CG"' > data/CpG/GSM1254259_HEK293-CT.CpG.txt

    # sum two strands for each CpG (so ignoring asymetric site)
    awk -v OFS='\t' '
      NR%2 {{ split($0, a) ; next }}
      {{print a[1], a[2], $2, a[7]+$7, a[8]+$8, (a[6]+$6)/2}}
    ' data/CpG/GSM1254259_HEK293-CT.CpG.txt > data/CpG/GSM1254259_HEK293-CT.CpG_Sum.txt

    # remove high copynumber cpgs
    awk '$6<1.5' data/CpG/GSM1254259_HEK293-CT.CpG_Sum.txt | cut -f 1-5 > data/CpG/GSM1254259_HEK293-CT.CpG_Sum_lowCN.txt

    # hg19 to hg38
    liftOver data/CpG/GSM1254259_HEK293-CT.CpG_Sum_lowCN.txt {input.chain} {output.cpg_hg38} data/CpG/GSM1254259_HEK293-CT.CpG.nonhg38.bed

    # Careful, liftover creates overlaping (duplicate CpG locations)

    # sort
    sort -k1,1 -k2,2n {output.cpg_hg38} > {output.cpg_hg38_sort}
    sort -k1,1 -k2,2n {input.windows} > {output.windows_namesort}

    # sum for each cpg Island
    bedtools map -c 4,5 -o sum -null 0 -a {input.cgi} -b {output.cpg_hg38_sort} > {output.cgimeth}

    # sum for each 100bp
    bedtools map -c 4,5 -o sum -null 0 -a {output.windows_namesort} -b {output.cpg_hg38_sort} > {output.100bp}

    awk -v OFS='\t' '{{print $1,$2,$3,($4+1)/($4+$5+2)}}' {output.cpg_hg38_sort} > {output.mean}

    # for IGV figure
    awk '$1=="chr1" && $2>53800000 && $2<54000000' {output.mean} > {output.mean_chr1}

    bedGraphToBigWig {output.mean_chr1} {input.chrom} {output.bigwig}

    # issue due to counting CpG as 1 position, but in meth counting it as 2bp long due to reverse strand, can have CpG==0 region with meth reads...
    #awk '$1=="chr1"  && $2>2442000' GSM1254259_HEK293-CT.CpG.hg38.sort.bed | head
    #cpg_Meth_100bp[chr=="chr1" & center_start==2442100]
    """

rule mappability:
  input:
    hg_sizes = "data/genome/hg38_sizes_23.chrom",
    zcwpw1_cin = "data/peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.bed"
  output:
    hg_sizes_alphasort = "data/genome/hg38_sizes_23_alphasort.chrom",
    bw = "data/mappability/k24.Umap.MultiTrackMappability.bw",
    bg = "data/mappability/k24.Umap.MultiTrackMappability.bedGraph",
    bed_keep = "data/mappability/24.Umap.MultiTrackMappability_Keep.bed",
    bed_exclude = "data/mappability/24.Umap.MultiTrackMappability_Exclude.bed",
    zcwpw1_cin_mappable = "data/mappability/Zcwpw1_peak_cin_24bpmappable.bed"
  shell:
    """
    # convert Hoffman lab mappability track to bed
    wget -P data/mappability/ https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/hg38/k24.Umap.MultiTrackMappability.bw

    ./bigWigToBedGraph {output.bw} {output.bg}

    # keep if mappabilty >0.75
    awk '$4 > 0.75' {bg} | LC_ALL=C sort -k1,1 -k2,2n -S5G --parallel=5 | \
    bedtools merge -i stdin > {output.bed_keep}

    sort -k1,1 {input.hg_sizes} > {output.hg_sizes_alphasort}

    # invert for exclusion bed
    awk '$1!="chrY"' {output.bed_keep} | \
      bedtools complement -i stdin -g {output.hg_sizes_alphasort} \
      > {output.bed_exclude}

    # censor Zcwpw1 peaks that aren't mappable if using 24bp (for overlap with Chip that did use 24bp..)
    bedtools slop -i {bed_exclude} -g {output.hg_sizes_alphasort} -b 10 | \
      bedtools subtract -A -a {input.zcwpw1_cin} -b stdin \
      > {output.zcwpw1_cin_mappable}
    """

rule download_ENCODE_beds:
  output:
    "data/ENCODE_beds/ENCFF{bedcode}.bed"
  shell:
    """
    wget -P data/ENCODE_beds/ "https://www.encodeproject.org/files/$i/@@download/ENCFF{output}.bed.gz"
    gz -k data/ENCODE_beds/ENCFF{output}.bed.gz
    """



ENCODE_BEDS = ["ENCFF314ZAL",
                "ENCFF338RAX",
                "ENCFF860DHS",
                "ENCFF959SPJ",
                "ENCFF131RPK",
                "ENCFF108BVL",
                "ENCFF418KWK",
                "ENCFF451UZW",
                "ENCFF464QPC",
                "ENCFF700RBU",
                "ENCFF483QXH",
                "ENCFF446OZF",
                "ENCFF786NME",
                "ENCFF478NGK",
                "ENCFF129ADK",
                "ENCFF204MYY",
                "ENCFF678TFE",
                "ENCFF081QMM",
                "ENCFF101AKQ",
                "ENCFF611CFB",
                "ENCFF342JBJ",
                "ENCFF538EDC",
                "ENCFF422AIH",
                "ENCFF295WQL",
                "ENCFF665UWA",
                "ENCFF439CWL",
                "ENCFF280SGN",
                "ENCFF653UGW",
                "ENCFF435BYC",
                "ENCFF567BLE",
                "ENCFF151LTX",
                "ENCFF403LQH",
                "ENCFF720TFF",
                "ENCFF108PRX"
                ]

rule analyse_peaks:
  input:
    "functions.R",
    expand("data/ENCODE_beds/{sample}.bed", sample=ENCODE_BEDS),
    rmd = "analysis/PeakCallingMA.Rmd",
    cpgi = "data/CpG/cpgIslandExtUnmasked.txt",
    repeats = "data/repeats/Repeat_MaskerSorted.bed",
    mappability = "../data/mappability/Zcwpw1_peak_cin_24bpmappable.bed"
  output:
    md = "results/PeakCallingMA.md"
  shell:
    """
    R -e "knitr::knit('{input.rmd}', '{output.md}')"
    """

rule analyse_100bp_windows:
  input:
    "peaks/ForceCalledPeaks_NA15-SRR5627150_AND_NA15-SRR5627151_vs_NA15-SRR5627142_AT_cpgIslandExtUnmasked.bed.bed",
    "peaks/ForceCalledPeaks_WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144_AT_cpgIslandExtUnmasked.bed.bed",
    "data/CpG/cpgIsland_Meth.bed",
    "data/CpG/CpG_hg38_100bpwindows.bed",
    "data/motifs/motif_locs_M7_100bpwindows.bed",
    "data/motifs/fimo_M4_WG_100bpwindows.bed",
    "data/motifs/motif_locs_CHB6_100bpwindows.bed"
    "data/peaks/Zcw_random_1bp.100bpWindows.bed",
    "data/pratto_dmc1/Pratto_human_DSB_Aintersect_autosomal_hg38.bed",
    "data/peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.100bpWindows.bed",
    "data/repeats/Repeat_MaskerSorted_NotAlu_100bp.bed",
    "data/repeats/Repeat_MaskerSorted_Alu_QC_100bp_f0.5.bed",
    "data/repeats/Repeat_MaskerSorted_Alu_QC_100bp.bed"
    "data/CpG/cpgIslandExtUnmasked_100bpwindows.bed",
    "data/CpG/cpgIslandExtUnmasked_100bpwindows_f1.bed",
    "data/CpG/cpg_Meth_100bpwindows.bed",
    "peaks/ForceCalledPeaks_NA15-SRR5627150_vs_NA15-SRR5627142_AT_genome.windows.100wide.100slide.bed.bed" # and others



