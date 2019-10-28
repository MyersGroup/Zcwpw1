############################
####### DMC1
############################

wget -p dmc1/ ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59836/suppl/GSE59836%5FPeak%5Fdata%5FSupplementary%5FFile%5F1%2Etxt%2Egz
gunzip dmc1/GSE59836_Peak_data_Supplementary_File_1.txt.gz

# subset to only AorB allele peaks observed in all samples & compute average strength
awk -v OFS='\t' '{ if ($16 == 1 && ($3-$2)<3000 && $1!="chrY" && $1!="chrX") { print $1, $2, $3, ($4+$5)/2, ($6+$7)/4 } }' \
  dmc1/GSE59836_Peak_data_Supplementary_File_1.txt \
  > dmc1/Pratto_human_DSB_Aintersect_autosomal_hg19.bed
# 18342
# 18326

liftOver dmc1/Pratto_human_DSB_Aintersect_autosomal_hg19.bed hg19ToHg38.over.chain.gz \
  dmc1/Pratto_human_DSB_Aintersect_autosomal_hg38.bed \
  dmc1/Pratto_human_DSB_Aintersect_autosomal_Nonhg38.bed

# awk -v OFS='\t' '{ if ($11 == 1 && ($3-$2)<3000 && $1!="chrY" && $1!="chrX") { print $1, $2, $3, ($4+$5)/2, ($6+$7)/4 } }' GSE59836_Peak_data_Supplementary_File_1.txt \ > Pratto_human_DSB_AB_autosomal_hg19.bed
# wc -l Pratto_human_DSB_AB_autosomal_hg19.bed
# #20336
#
# awk -v OFS='\t' '{ if ($9 == 1 && $10 == 1 && ($3-$2)<3000 && $1!="chrY" && $1!="chrX") { print $1, $2, $3, ($4+$5)/2, ($6+$7)/4 } }' GSE59836_Peak_data_Supplementary_File_1.txt \ > Pratto_human_DSB_AA_autosomal_hg19.bed
# #21699


############################
####### Repeats
############################

curl 'http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=738896965_zTgYMjAIHdEHcvRZcCteap7Zqbpt&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_rmsk&hgta_ctDesc=table+browser+query+on+rmsk&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED' \
  -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.14; rv:67.0) Gecko/20100101 Firefox/67.0' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
  -H 'Accept-Language: en-GB,en;q=0.5' \
  --compressed -H 'Connection: keep-alive' \
  -H 'Referer: http://genome.ucsc.edu/cgi-bin/hgTables' \
  -H 'Cookie: hguid.genome-euro=362527140_Mi0jIpneXLkgeplCtz5aa6JO13ha; hguid=673678055_VvqPPBb5iFgOJEF3WMG7M3t7iual' \
  -H 'Upgrade-Insecure-Requests: 1' > repeats/Repeat_Masker.bed.gz

gunzip -d repeats/Repeat_Masker.bed.gz

grep -P 'chr[0-9XY]+\t' repeats/Repeat_Masker.bed | sort -k1,1 -k2,2n  > repeats/Repeat_MaskerSorted.bed

# split into Alu and not
grep 'Alu[a-zA-Z]' repeats/Repeat_MaskerSorted.bed > repeats/Repeat_MaskerSorted_Alu.bed
grep -v 'Alu[a-zA-Z]' repeats/Repeat_MaskerSorted.bed > repeats/Repeat_MaskerSorted_NotAlu.bed

# check Alu with no family name not present
grep -P 'Alu\t' repeats/Repeat_MaskerSorted_Alu.bed | head

awk '$3-$2 > 250 && $3-$2 < 350' repeats/Repeat_MaskerSorted_Alu.bed > repeats/Repeat_MaskerSorted_Alu_QC.bed


bedtools getfasta -s -fi motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed repeats/Repeat_MaskerSorted_Alu_QC.bed -name > repeats/AluRepeats.fa


############################
####### Check Enrichment at Alus
############################

awk '$3-$2 == 299' repeats/Repeat_Masker_Alu.bed  > Repeat_Masker_Alu_299.bed

bwtool aggregate 300:300:300 Repeat_Masker_Alu_299.bed bedgraphs/depth_WTCHG_538916_221156.bigWig Alu_meta299.txt -fill=0
bwtool aggregate 300:300:300 Repeat_Masker_Alu_299.bed bedgraphs/depth_WTCHG_538916_217108.bigWig Alu_meta_In299.txt -fill=0

bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_NA15-SRR5627148.bigWig repeats/Alu_UT_H3K36me3.txt -fill=0
bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_NA15-SRR5627150.bigWig repeats/Alu_UT_H3K4me3.txt -fill=0
bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_NA15-SRR5627142.bigWig repeats/Alu_UT_input.txt -fill=0




############################
####### Isoforms from GTEX data
############################
cd GTEX
wget -P GTEX/ https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
zgrep -E "ENSG00000078487.17|transcript" GTEX/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz > GTEX/zcwpw1_isoforms.tsv

wget -P GTEX/ https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt



#####################
######## CpGs #######
#####################

mkdir CpG

# Random Zcw Peaks
head -23 ../single-cell/sequencing/metadata/hg38_sizes.chrom > hg38_sizes23.bed
bedtools random -n 1264519 -l 300 -g hg38_sizes23.bed -seed 71346 | sort -k1,1 -k2,2n > peaks/Zcw_random.bed

# # all CpG in Chr1
# head -1 ../single-cell/sequencing/metadata/hg38_sizes.bed |
# bedtools getfasta -fi motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed stdin |
# ./seqkit locate --ignore-case --non-greedy --only-positive-strand -p CG --bed |
# sed s/:0-248956422//g > CpG_chr1.bed

# All CpG positions in genome
./seqkit locate --ignore-case --non-greedy --only-positive-strand -p CG --bed motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa > CpG/CpG_hg38.bed

cut -f 1-3 CpG/CpG_hg38.bed > CpG/CpG_hg38_cut.bed

# count CpG at each peak (Chr1)
# tail -n +2 peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.bed |
# bedtools slop -i stdin -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 150 | grep -P "chr1\t" |
# bedtools intersect -a stdin -b CpG_chr1.bed -c > peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL_wCpGcount_chr1.bed

# count CpG at each peak Whole genome
tail -n +2 peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.bed |
bedtools slop -i stdin -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 150 |
bedtools intersect -a stdin -b CpG/CpG_hg38.bed -c > peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL_wCpGcount.bed

# and for CVC peaks
tail -n +2 peaks/SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL.bed |
bedtools slop -i stdin -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 150 |
bedtools intersect -a stdin -b CpG/CpG_hg38.bed -c > peaks/SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL_wCpGcount.bed


# Count CpG at random Zcw peaks
bedtools intersect -a peaks/Zcw_random.bed -b CpG/CpG_hg38.bed -c > peaks/Zcw_random_CpGcount.bed

# 115320 peaks, each 300bp
# 2316765 CpG positions
# chr1 248956422
# 2316765 / (248956422/300) = 2.79 CpG per 300bp
# but clustering
# try randomisation instead
# head -1 ../single-cell/sequencing/metadata/hg38_sizes.chrom > chr1_size.bed
# bedtools random -n 115320 -l 300 -g chr1_size.bed -seed 71346 | sort -k1,1 -k2,2n > chr1_random.bed
#
# bedtools intersect -a chr1_random.bed -b CpG_chr1.bed -c > chr1_random_CpGcount.bed


####################################
###### Count CpG & CpG Islands per 100bp window
####################################

bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b CpG/CpG_hg38_cut.bed -c > CpG/CpG_hg38_100bpwindows.bed

##### CpG Islands per 100bp

wget -P CpG/ http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz
gzip -d CpG/cpgIslandExtUnmasked.txt.gz
cut -f 2-4 CpG/cpgIslandExtUnmasked.txt > CpG/cpgIslandExtUnmasked.bed

bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b CpG/cpgIslandExtUnmasked.bed -f 0.1 -c > CpG/cpgIslandExtUnmasked_100bpwindows.bed

bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b CpG/cpgIslandExtUnmasked.bed -f 1 -c > CpG/cpgIslandExtUnmasked_100bpwindows_f1.bed


####################################
###### Bisulphite sequencing data
####################################

#https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000533
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1254259
wget -P CpG/ ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1254nnn/GSM1254259/suppl/GSM1254259%5FHEK293%2DCT%2Ecout%2Etxt%2Egz

# only CpG meth
zcat CpG/GSM1254259_HEK293-CT.cout.txt.gz | awk ' $4 == "CG"' > CpG/GSM1254259_HEK293-CT.CpG.txt

# sum two strands for each CpG (so ignoring asymetric site)
awk -v OFS='\t' '
  NR%2 { split($0, a) ; next }
  {print a[1], a[2], $2, a[7]+$7, a[8]+$8, (a[6]+$6)/2}
' CpG/GSM1254259_HEK293-CT.CpG.txt > CpG/GSM1254259_HEK293-CT.CpG_Sum.txt

# remove high copynumber cpgs
awk '$6<1.5' CpG/GSM1254259_HEK293-CT.CpG_Sum.txt | cut -f 1-5 > CpG/GSM1254259_HEK293-CT.CpG_Sum_lowCN.txt

# hg19 to hg38
liftOver CpG/GSM1254259_HEK293-CT.CpG_Sum_lowCN.txt hg19ToHg38.over.chain.gz CpG/GSM1254259_HEK293-CT.CpG.hg38.bed CpG/GSM1254259_HEK293-CT.CpG.nonhg38.bed

# Careful, liftover creates overlaping (duplicate CpG locations)

# sort
sort -k1,1 -k2,2n CpG/GSM1254259_HEK293-CT.CpG.hg38.bed > CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.bed
sort -k1,1 -k2,2n forcepeaks/genome.windows.100wide.100slide.bed > forcepeaks/genome.windows.100wide.100slide.namesort.bed

# sum for each cpg Island
bedtools map -c 4,5 -o sum -null 0 -a CpG/cpgIslandExtUnmasked.bed -b CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.bed > CpG/cpgIsland_Meth.bed

# sum for each 100bp
bedtools map -c 4,5 -o sum -null 0 -a forcepeaks/genome.windows.100wide.100slide.namesort.bed -b CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.bed > CpG/cpg_Meth_100bpwindows.bed

awk -v OFS='\t' '{print $1,$2,$3,($4+1)/($4+$5+2)}' CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.bed > CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.bed

# for IGV figure
awk '$1=="chr1" && $2>53800000 && $2<54000000' CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.bed > CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.chr1.bed

bedGraphToBigWig CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.chr1.bed ../single-cell/sequencing/metadata/hg38_sizes.chrom CpG/GSM1254259_HEK293-CT.CpG.hg38.sort.mean.chr1.bigWig

# issue due to counting CpG as 1 position, but in meth counting it as 2bp long due to reverse strand, can have CpG==0 region with meth reads...
#awk '$1=="chr1"  && $2>2442000' GSM1254259_HEK293-CT.CpG.hg38.sort.bed | head
#cpg_Meth_100bp[chr=="chr1" & center_start==2442100]


############################
####### Count per 100bp
############################

# Count Zcw peaks per 100bp (max 1 as min sep=250)
bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.bed -c > peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.100bpWindows.bed

# random peaks (made above) per 100bp
awk -v OFS='\t' '{print $1, $2, $2+1}' peaks/Zcw_random.bed > peaks/Zcw_random_1bp.bed

bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b peaks/Zcw_random_1bp.bed -c > peaks/Zcw_random_1bp.100bpWindows.bed


# Count Alus per 100bp windows

# f= fraction of A that must overlap

bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b repeats/Repeat_MaskerSorted_Alu_QC.bed -f 0.1 -c > repeats/Repeat_MaskerSorted_Alu_QC_100bp.bed

bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b repeats/Repeat_MaskerSorted_Alu_QC.bed -f 1 -c > repeats/Repeat_MaskerSorted_Alu_QC_100bp_f0.5.bed

bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b repeats/Repeat_MaskerSorted_NotAlu.bed -f 0.1 -c > repeats/Repeat_MaskerSorted_NotAlu_100bp.bed




############################
####### ENCODE beds
############################



declare -a arr=("ENCFF314ZAL"
                "ENCFF338RAX"
                "ENCFF860DHS"
                "ENCFF959SPJ"
                "ENCFF131RPK"
                "ENCFF108BVL"
                "ENCFF418KWK"
                "ENCFF451UZW"
                "ENCFF464QPC"
                "ENCFF700RBU"
                "ENCFF483QXH"
                "ENCFF446OZF"
                "ENCFF786NME"
                "ENCFF478NGK"
                "ENCFF129ADK"
                "ENCFF204MYY"
                "ENCFF678TFE"
                "ENCFF081QMM"
                "ENCFF101AKQ"
                "ENCFF611CFB"
                "ENCFF342JBJ"
                "ENCFF538EDC"
                "ENCFF422AIH"

                "ENCFF295WQL"
                "ENCFF665UWA"
                "ENCFF439CWL"
                "ENCFF280SGN"
                "ENCFF653UGW"
                "ENCFF435BYC"
                "ENCFF567BLE"
                "ENCFF151LTX"
                "ENCFF403LQH"
                "ENCFF720TFF"
                "ENCFF108PRX"
                )

## now loop through the above array
for i in "${arr[@]}"
do
    wget -P ENCODE_beds/ "https://www.encodeproject.org/files/$i/@@download/$i.bed.gz"
done


############################
####### SPO11
############################


sed 's/chr//' B6_Spo11.bedgraph > B6_Spo11_clean.bedgraph

zgrep -v '_' B6_Spo11.bedgraph | grep -vP 'M|X|Y' | sed 's/chr//' > B6_Spo11_clean.bedgraph

bedGraphToBigWig B6_Spo11_clean.bedgraph ../../single-cell/sequencing/metadata/mm10_sizes.chrom B6_Spo11.bedgraph.bigWig

bwtool matrix -fill=0 -decimals=1 -tiled-averages=5 5000:5000 B6.bed B6_Spo11.bedgraph.bigWig B6_Spo11_atB6.bwm


############################
####### Mappability
############################


# Mappability track from genome browser (>7GB)
# https://epgg-test.wustl.edu/d/hg38/hg38.mappability.75.bigwig

# convert Hoffman lab mappability track to bed
wget -P mappability/ https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/hg38/k24.Umap.MultiTrackMappability.bw

./bigWigToBedGraph mappability/k24.Umap.MultiTrackMappability.bw mappability/k24.Umap.MultiTrackMappability.bedGraph

# keep if mappabilty >0.75
awk '$4 > 0.75' mappability/k24.Umap.MultiTrackMappability.bedGraph | LC_ALL=C sort -k1,1 -k2,2n -S5G --parallel=5 | \
  bedtools merge -i stdin > mappability/24.Umap.MultiTrackMappability_Keep.bed

sort -k1,1 ../hg38_sizes_23.chrom > ../hg38_sizes_23_alphasort.chrom

# invert for exclusion bed
awk '$1!="chrY"' mappability/24.Umap.MultiTrackMappability_Keep.bed | \
  bedtools complement -i stdin -g ../hg38_sizes_23_alphasort.chrom \
  > mappability/24.Umap.MultiTrackMappability_Exclude.bed

# censor Zcwpw1 peaks that aren't mappable if using 24bp (for overlap with Chip that did use 24bp..)
bedtools slop -i mappability/24.Umap.MultiTrackMappability_Exclude.bed -g hg38_sizes_23_alphasort.chrom -b 10 | \
  bedtools subtract -A -a peaks/SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.bed -b stdin \
  > mappability/Zcwpw1_peak_cin_24bpmappable.bed


bwtool extract bed ch6_sample.bed mappability/k24.Umap.MultiTrackMappability.bw mappability/k24.Umap.MultiTrackMappability_bwtool.bed

bedtools merge -d 10 -i mappability/24.Umap.MultiTrackMappability.bed > mappability/24.Umap.MultiTrackMappability_d10.bed

#

wget -P mappability/ https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/hg38/k36.Umap.MultiTrackMappability.bw
./bigWigToBedGraph mappability/k36.Umap.MultiTrackMappability.bw mappability/k36.Umap.MultiTrackMappability.bedGraph
awk '$4 > 0.8' mappability/k36.Umap.MultiTrackMappability.bedGraph | LC_ALL=C sort -k1,1 -k2,2n -S5G --parallel=5 | bedtools merge -i stdin > mappability/36.Umap.MultiTrackMappability_Keep.bed
bedtools merge -d 10 -i mappability/36.Umap.MultiTrackMappability.bed > mappability/36.Umap.MultiTrackMappability_d10.bed

#

# example - not sure it's correct due to large streaches
wget -P mappability/ https://github.com/xuefzhao/Reference.Mappability/raw/master/hg38.50mer/hg38.50mer.Excludable.bed


sed 's/ /\t/g' mappability/hg38.50mer.Excludable.bed |
bedtools intersect -v -a  peaks/SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL.bed \
                          -b stdin > mappability/SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_217108_AND_WTCHG_538916_220144.p0.000001.sep250.ALL_QCMappability50.bed


# Use GEM
# Create mappability track with GEM
#(only exists publicly on hg19)

wget -P mappability/ https://downloads.sourceforge.net/project/gemlibrary/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2

tar -jxvf mappability/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2

export PATH=/homes/wells/saxony/zcwpw1/mappability/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/:$PATH

wget https://github.com/xuefzhao/Reference.Mappability/raw/master/Scripts/bedGraphTobed
chmod 777 bedGraphTobed
bedGraphTobed

# must change from py2 to py3 add paren to print statements, and file arg instead of >>


# unzip fasta reference
zcat /homes/wells/saxony/single-cell/sequencing/metadata/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz > /homes/wells/saxony/single-cell/sequencing/metadata/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

# create index
gem-indexer -i /homes/wells/saxony/single-cell/sequencing/metadata/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna  -o /homes/wells/saxony/single-cell/sequencing/metadata/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_setGEM

# generate mappability scores
gem-mappability -I /homes/wells/saxony/single-cell/sequencing/metadata/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_setGEM.gem -l 75 -o mappability/hg38.75mer -T 12
gem-2-wig -I /homes/wells/saxony/single-cell/sequencing/metadata/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_setGEM.gem -i mappability/hg38.75mer.mappability -o mappability/hg38.75mer

# remove erreneous names in chromosome names
sed "s/ AC//" mappability/hg38.75mer.wig > mappability/hg38.75mer2.wig

# convert to bed
./wigToBigWig mappability/hg38.75mer2.wig /homes/wells/saxony/single-cell/sequencing/metadata/hg38_sizes.chrom mappability/hg38.75mer.bw
./bigWigToBedGraph mappability/hg38.75mer.bw mappability/hg38.75mer.bedGraph
./bedGraphTobed mappability/hg38.75mer.bedGraph mappability/hg38.75mer.Includable.bed 0.33
../bedGraphTobed mappability/hg38.75mer.bedGraph mappability/hg38.75mer.Includable.t0.5.bed 0.5

##### NOT EXCLUDABLE BUT INCLUDABLE ######
#find -type f -name '*' -execdir rename -v 's:Excludable:Includable:g' {} \;

# subset smaller files
grep -P 'chr6 ' mappability/hg38.75mer.Includable.bed | sed -e 's/ /\t/g' > mappability/hg38.75mer.Includable.chr6.bed
grep -P 'chr1 ' mappability/hg38.75mer.Includable.bed | sed -e 's/ /\t/g' > mappability/hg38.75mer.Includable.chr1.bed
grep -P 'chr1 ' mappability/hg38.75mer.Includable.t0.5.bed | sed -e 's/ /\t/g' > mappability/hg38.75mer.Includable.t0.5.chr1.bed

head -1 ../single-cell/sequencing/metadata/hg38_sizes.chrom > mappability/hg38_chr1.chrom

bedtools complement -i mappability/hg38.75mer.Includable.chr1.bed -g ../hg38_chr1.chrom > mappability/hg38.75mer.Excludable.chr1.bed
bedtools complement -i mappability/hg38.75mer.Includable.t0.5.chr1.bed -g ../hg38_chr1.chrom > mappability/hg38.75mer.Excludable.t0.5.chr1.bed

bedtools merge -i mappability/hg38.75mer.Includable.chr1.bed -d 100 > mappability/hg38.75mer.Includable.chr1.merge100.bed
bedtools merge -i mappability/hg38.75mer.Excludable.chr1.bed -d 500 > mappability/hg38.75mer.Excludable.chr1.merge500.bed
bedtools merge -i mappability/hg38.75mer.Excludable.t0.5.chr1.bed -d 500 > mappability/hg38.75mer.Excludable.t0.5.chr1.merge500.bed

rg 'chr6\t' mappability/hg38.75mer.bedGraph > mappability/hg38.75mer.chr6.bedGraph
rg 'chr1\t' mappability/hg38.75mer.bedGraph > mappability/hg38.75mer.chr1.bedGraph


############################
####### Other
############################


# small bam file for sample
#chr6:107,211,482-107,264,849
echo "chr6  107211482  107264849" ch6_sample.bed
samtools view filtered/WTCHG_538916_221156.bam -L ch6_sample.bed -b -o ch6_sample.bam
samtools index ch6_sample.bam



# Hells

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM835nnn/GSM835828/suppl/GSM835828%5Fwt%5Fhells%2Ebed%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM835nnn/GSM835828/suppl/GSM835828%5Fwt%5Fhells%5Fpeaks%2Ebed%2Egz

wget https://raw.githubusercontent.com/arq5x/bedtools/master/genomes/mouse.mm9.genome
zcat GSM835828_wt_hells.bed.gz | sort -k 1,1 > GSM835828_wt_hells.sorted.bed
bedtools genomecov -bg -i GSM835828_wt_hells.sorted.bed -g mouse.mm9.genome > GSM835828_wt_hells.bedGraph

scp wells@stats2:saxony/zcwpw1/GSM835828_wt_hells.bedGraph .



# Remove more alignments from BWA MEM tags
```{r}
sambamba view -t 12 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" filtered/WTCHG_538916_221156.bam -o WTCHG_538916_221156_SBBF.bam

samtools flagstat -@ 5 WTCHG_538916_221156_SBBF.bam
# 80705408 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 80705408 + 0 mapped (100.00% : N/A)
# 80705408 + 0 paired in sequencing
# 40448867 + 0 read1
# 40256541 + 0 read2
# 80705408 + 0 properly paired (100.00% : N/A)
# 80705408 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)


samtools flagstat -@ 5 filtered/WTCHG_538916_221156.bam
# 88007214 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 88007214 + 0 mapped (100.00% : N/A)
# 88007214 + 0 paired in sequencing
# 44003607 + 0 read1
# 44003607 + 0 read2
# 88007214 + 0 properly paired (100.00% : N/A)
# 88007214 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
```

# test MDB4 & SETDB1
```{r}

conda install -c bioconda ucsc-bigwigtowig
https://github.com/ENCODE-DCC/kentUtils/blob/v302.1.0/bin/linux.x86_64/bigWigToWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig

wget https://www.encodeproject.org/files/ENCFF000ZTL/@@download/ENCFF000ZTL.bigWig
./bigWigToWig -chrom=chr7 -start=139000000 -end=139600000 ENCFF000ZTL.bigWig ENCFF000ZTL.Wig
liftOver ENCFF000ZTL.Wig hg19ToHg38.over.chain.gz ENCFF000ZTL_hg38.Wig ENCFF000ZTL_Nonhg38.Wig
./wigToBigWig ENCFF000ZTL_hg38.Wig ../single-cell/sequencing/metadata/hg38_sizes.chrom ENCFF000ZTL_hg38.BigWig


wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1010nnn/GSM1010740/suppl/GSM1010740%5Fhg19%5FwgEncodeHaibTfbsHepg2Mbd4sc271530V0422111PkRep1%2EbroadPeak%2Egz

zcat GSM1010740_hg19_wgEncodeHaibTfbsHepg2Mbd4sc271530V0422111PkRep1.broadPeak.gz | cut -f 1-3 > MBD4.bed
liftOver MBD4.bed hg19ToHg38.over.chain.gz MBD4_hg38.bed MBD4_nothg38.bed

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1010nnn/GSM1010740/suppl/GSM1010740%5Fhg19%5FwgEncodeHaibTfbsHepg2Mbd4sc271530V0422111RawRep1%2EbigWig
./bigWigToWig -chrom=chr7 -start=139000000 -end=139600000 GSM1010740_hg19_wgEncodeHaibTfbsHepg2Mbd4sc271530V0422111RawRep1.bigWig MDB4.Wig
liftOver MDB4.Wig hg19ToHg38.over.chain.gz MDB4_hg38.Wig MDB4_nonhg38.Wig
./wigToBigWig MDB4_hg38.Wig ../single-cell/sequencing/metadata/hg38_sizes.chrom MDB4_hg38.BigWig
```


##### DeNovo Motif

# count seqkit located motifs by 100bp
bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b motifs/motif_locs_M7.bed -c > motifs/motif_locs_M7_100bpwindows.bed
bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b motifs/motif_locs_CHB6.bed -c > motifs/motif_locs_CHB6_100bpwindows.bed

# count fimo located motifs by 100bp
awk '$5>13' motifs/fimo_M4_WG.tsv | cut -f 1-3 | tail -n +2 |
bedtools intersect -a forcepeaks/genome.windows.100wide.100slide.bed -b stdin -c > motifs/fimo_M4_WG_100bpwindows.bed



