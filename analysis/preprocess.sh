
############################
####### Check Enrichment at Alus
############################

awk '$3-$2 == 299' repeats/Repeat_Masker_Alu.bed  > repeats/Repeat_Masker_Alu_299.bed

bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_WTCHG_538916_221156.bigWig repeats/Alu_meta299.txt -fill=0
bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_WTCHG_538916_217108.bigWig repeats/Alu_meta_In299.txt -fill=0

bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_NA15-SRR5627148.bigWig repeats/Alu_UT_H3K36me3.txt -fill=0
bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_NA15-SRR5627150.bigWig repeats/Alu_UT_H3K4me3.txt -fill=0
bwtool aggregate 300:300:300 repeats/Repeat_Masker_Alu_299.bed bedgraphs/depth_NA15-SRR5627142.bigWig repeats/Alu_UT_input.txt -fill=0


############################
####### Mappability
############################


# Mappability track from genome browser (>7GB)
# https://epgg-test.wustl.edu/d/hg38/hg38.mappability.75.bigwig




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



