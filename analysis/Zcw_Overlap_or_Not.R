
# split HP9N into overlap with ZcwCVC and not overlap

cut -f 1,4,5 motifs/SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL_QCfiltered.bed |
  bedtools slop -i stdin -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 100 |
  bedtools intersect -v -a deeptools/beds/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.bed \
-b stdin > HP9N_NOToverlap_ZcwCVC.bed

cut -f 1,4,5 motifs/SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL_QCfiltered.bed |
  bedtools slop -i stdin -g ../single-cell/sequencing/metadata/hg38_sizes.chrom -b 100 |
  bedtools intersect -u -a deeptools/beds/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.bed \
-b stdin > HP9N_overlap_ZcwCVC.bed


# sort by H3K4me3 on one side

computeMatrix reference-point \
-S deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
-R HP9N_overlap_ZcwCVC.bed \
-b 1000 \
-p 15 \
--skipZeros \
-o deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_sort_overlapZcw.gz \
--sortRegions descend \
--outFileSortedRegions deeptools/beds/HP9N_overlap_ZcwCVC_sorted.bed

## extract regions data

computeMatrix reference-point \
-S deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
deeptools/bigwigs/SRA_Altemose2015_SRR5627149_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
-R deeptools/beds/HP9N_overlap_ZcwCVC_sorted.bed \
-a 2000 \
-b 2000 \
-p 15 \
--skipZeros \
--sortRegions keep \
-o deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_overlapZcw.gz

## plot

plotHeatmap -m deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_overlapZcw.gz \
-out deeptools/plots/Methylation_AT_HP9N_overlapZcw.png \
--colorMap RdBu \
--samplesLabel H3K4me3, H3K36me3 \
--zMax 0.65 0.15 \
--yMax 0.45 0.1 \
--sortRegions keep \
-z regions \
--dpi 600 \
--outFileNameMatrix deeptools/plots/Methylation_AT_HP9N_overlapZcw.tab

# Not overlaping

# sort by H3K4me3 on one side

computeMatrix reference-point \
-S deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
-R HP9N_NOToverlap_ZcwCVC.bed \
-b 1000 \
-p 15 \
--skipZeros \
-o deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_sort_NOToverlapZcw.gz \
--sortRegions descend \
--outFileSortedRegions deeptools/beds/HP9N_NOToverlap_ZcwCVC_sorted.bed

## extract regions data

computeMatrix reference-point \
-S deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
deeptools/bigwigs/SRA_Altemose2015_SRR5627149_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
-R deeptools/beds/HP9N_NOToverlap_ZcwCVC_sorted.bed \
-a 2000 \
-b 2000 \
-p 15 \
--skipZeros \
--sortRegions keep \
-o deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_NOToverlapZcw.gz

## plot

plotHeatmap -m deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_NOToverlapZcw.gz \
-out deeptools/plots/Methylation_AT_HP9N_NOToverlapZcw.png \
--colorMap RdBu \
--samplesLabel H3K4me3, H3K36me3 \
--zMax 0.25 0.1 \
--yMax 0.15 0.05 \
--sortRegions keep \
-z regions \
--dpi 600 \
--sortRegions keep \
--outFileNameMatrix deeptools/plots/Methylation_AT_HP9N_NOToverlapZcw.tab



# Profile plots only, overlaping vs not

computeMatrix reference-point \
-S deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
deeptools/bigwigs/SRA_Altemose2015_SRR5627149_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
-R HP9N_NOToverlap_ZcwCVC.bed \
HP9N_overlap_ZcwCVC.bed \
-a 2000 \
-b 2000 \
-p 15 \
--skipZeros \
-o deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_overlapZcw_andNOT.gz

plotProfile    \
-m deeptools/matrices/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Heatmap_overlapZcw_andNOT.gz \
-out deeptools/plots/overlapvsNot.pdf \
--refPointLabel "Motif Center" \
--regionsLabel OverlapZcw NotOverlapZcw \
--samplesLabel H3K4 H3K36 \
--numPlotsPerRow 4 \
--yMax 0.4 0.1
