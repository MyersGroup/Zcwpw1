

# awk -v OFS='\t' '{print $1,$2,$3,"0","0",$6; }' HP9N_overlap_ZcwCVC.bed > HP9N_overlap_ZcwCVC.bed6
# awk -v OFS='\t' '{print $1,$2,$3,"0","0",$6; }' HP9N_NOToverlap_ZcwCVC.bed > HP9N_NOToverlap_ZcwCVC.bed6
#
# awk -v OFS='\t' '{print $1,$2,$3,"0","0",$6; }' deeptools/beds/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.bed  > HP9N_MotifCenterStranded.bed6

# bwtool matrix 2000:2000 HP9N_overlap_ZcwCVC.bed6 bedgraphs/depth_WTCHG_538916_223180.bigWig,bedgraphs/depth_WTCHG_538916_221156.bigWig bwmatrices/depth_WTCHG_538916_223180_HP9N_overlap_ZcwCVC.bwm -fill=0
# bwtool matrix -fill=0 -decimals=0 2000:2000 HP9N_overlap_ZcwCVC.bed6 bedgraphs/depth_WTCHG_538916_223180.bigWig,bedgraphs/depth_WTCHG_538916_221156.bigWig bwmatrices/depth_WTCHG_538916_223180_HP9N_overlap_ZcwCVC.bwm
# bwtool matrix -fill=0 -decimals=1 -tiled-averages=5 2000:2000 HP9N_overlap_ZcwCVC.bed6 bedgraphs/depth_WTCHG_538916_223180.bigWig,bedgraphs/depth_WTCHG_538916_221156.bigWig bwmatrices/depth_WTCHG_538916_223180_HP9N_overlap_ZcwCVC.bwm



# SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.bed00 # Nterm
# SingleBasePeaks.SRA_Altemose2015_SRR5627146_AND_SRA_Altemose2015_SRR5627147_vs_SRA_Altemose2015_SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded.bed00 # C term combo

# Zcw_HP9C & Zcw
#bwtool matrix -fill=0 -decimals=1 -tiled-averages=5 2000:2000 bed6/SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded_Q00.bw.bed00 bedgraphs/depth_WTCHG_538916_223180.bigWig bwmatrices/depth_WTCHG_538916_223180_HP9N.bwm

snakemake --cores 15 --snakefile pipelines/Plot_Heatmap.py -npr

library(data.table)
library(ComplexHeatmap)


normalise_matrix <- function(sample, input, pseudocount=1, range=c(0.025, 1-0.025)){
  norm <- (sample+pseudocount)/(input+pseudocount)
  # 0/0 = 0
  norm[is.na(norm)] <- 0

  # x/0 = Inf, pass
  norm[!is.finite(norm)] <- NA

  norm[norm>quantile(norm, na.rm=T, range[2])] <- quantile(norm, na.rm=T, range[2])
  norm[norm<quantile(norm, na.rm=T, range[1])] <- quantile(norm, na.rm=T, range[1])

  return(norm)
}




makeColumn <- function(sample, input, title, ordering){
  profile <- colMeans(sample, na.rm = T) / colMeans(input, na.rm = T)
  profileplot = HeatmapAnnotation("Avg"=anno_lines(profile,
                                                   height = unit(4, "cm"),
                                                   axis_param = list(facing="inside"),
                                                   ylim = c(min(1, min(profile)*0.9),
                                                            max(2,max(profile)*1.1))),
                                  show_annotation_name=FALSE)

  normalised <- normalise_matrix(sample,input)

  heatmapcolumn =  ComplexHeatmap::Heatmap(log(normalised)[ordering,],
                                           cluster_rows = F,
                                           cluster_columns = F,
                                           name = title,
                                           column_title = title,
                                           top_annotation = profileplot,
                                           show_column_names = F,
                                           col=RColorBrewer::brewer.pal(11,"RdBu"),
                                           column_split = rep(c("L","R"), each = ncol(normalised)/2),
                                           column_gap = unit(0, "mm"),
                                           border = TRUE,
                                           column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                           use_raster = TRUE,
                                           raster_quality = 8,
                                           heatmap_legend_param = list(direction = "horizontal"))

  return(heatmapcolumn)
}

plotHeatmap <- function(locations, name){

  Zcw <- as.matrix(fread(paste0("bwmatrices/WTCHG_538916_221156_AT_",locations,".bwm")))
  Zcw_hP9C <- as.matrix(fread(paste0("bwmatrices/WTCHG_538916_223180_AT_",locations,".bwm")))
  Zcw_In <- as.matrix(fread(paste0("bwmatrices/WTCHG_538916_220144_AT_",locations,".bwm")))
  Zcw_hP9C_In <- as.matrix(fread(paste0("bwmatrices/WTCHG_538916_217108_AT_",locations,".bwm")))

  H3K4_hP9C <- as.matrix(fread(paste0("bwmatrices/NA15-SRR5627152_AND_NA15-SRR5627153_AT_",locations,".bwm")))
  H3K36_hP9C <- as.matrix(fread(paste0("bwmatrices/NA15-SRR5627149_AT_",locations,".bwm")))
  Input_hP9C <- as.matrix(fread(paste0("bwmatrices/NA15-SRR5627143_AT_",locations,".bwm")))

  hP9C <- as.matrix(fread(paste0("bwmatrices/NA15-SRR5627146_AT_",locations,".bwm")))

  UntIn <- as.matrix(fread(paste0("bwmatrices/NA15-SRR5627142_AT_",locations,".bwm")))
  UntH3K4 <- as.matrix(fread(paste0("bwmatrices/NA15-SRR5627150_AND_NA15-SRR5627151_AT_",locations,".bwm")))
  UntH3K36 <- as.matrix(fread(paste0("bwmatrices/NA15-SRR5627148_AT_",locations,".bwm")))

  # ordering <- order(apply(H3K4_hP9C_vIn, 1, which.max)) # induces artificial ridge

  hP9C_vIn <- normalise_matrix(hP9C,Input_hP9C)
  ordering <- order(-apply(hP9C_vIn[,300:500], 1, function(x) mean(x, na.rm =T)))

  ht_list <- makeColumn(Zcw, Zcw_In+Zcw_hP9C_In, "ZHA vs In", ordering) +
    makeColumn(Zcw_hP9C, Zcw_In+Zcw_hP9C_In, "ZHA_hP9V5 vs In", ordering) +
    makeColumn(Zcw_hP9C, Zcw, "ZHA_hP9V5 vs ZHA", ordering) +
    makeColumn(H3K4_hP9C, Input_hP9C, "H3K4_hP9HA vs InhP9HA", ordering) +
    makeColumn(H3K36_hP9C, Input_hP9C, "H3K36_hP9HA vs InhP9HA", ordering) +
    makeColumn(hP9C, Input_hP9C, "hP9HA vs InhP9HA", ordering) +
    makeColumn(UntH3K4,UntIn, "UntH3K4 vs InhUnt", ordering) +
    makeColumn(UntH3K36,UntIn, "UntH3K36 vs InUnt", ordering)

  pdf(width = 15, height = 10, paste0("bwplots/NormalisedhP9C_vInOrder_",name,"-locations.pdf"))
  draw(ht_list, merge_legend=T, heatmap_legend_side="bottom")
  dev.off()

}

plotHeatmap("SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded_Q00","HP9C")
plotHeatmap("SingleBasePeaks.NA15-SRR5627138_AND_NA15-SRR5627139_vs_NA15-SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Q00","HP9N")
plotHeatmap("SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL_Q00","CHP9C")
plotHeatmap("SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL_Q00","hZcwCTvsST")
plotHeatmap("SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108.p0.000001.sep250.ALL_Q00","hZcwSTvsIn")



#### Order by left right H3K4

ordering <- order(apply(H3K4_hP9C[,1:400], 1, function(x) mean(x, na.rm =T)) / apply(H3K4_hP9C[,401:800], 1, function(x) mean(x, na.rm =T)))

ht_list <- makeColumn(Zcw, Zcw_In+Zcw_hP9C_In, "ZHA vs In", ordering) +
  makeColumn(Zcw_hP9C, Zcw_In+Zcw_hP9C_In, "ZHA_hP9V5 vs In", ordering) +
  makeColumn(Zcw_hP9C, Zcw, "ZHA_hP9V5 vs ZHA", ordering) +
  makeColumn(H3K4_hP9C, Input_hP9C, "H3K4_hP9HA vs InhP9HA", ordering) +
  makeColumn(H3K4b_hP9C, Input_hP9C, "H3K4_hP9HA Rep2 vs InhP9HA", ordering) +
  makeColumn(H3K36_hP9C, Input_hP9C, "H3K36_hP9HA vs InhP9HA", ordering) +
  makeColumn(hP9C, Input_hP9C, "hP9HA vs InhP9HA", ordering) +
  makeColumn(UntH3K4,UntIn, "UntH3K4 vs InhUnt", ordering) +
  makeColumn(UntH3K36,UntIn, "UntH3K36 vs InUnt", ordering)

pdf(width = 15, height = 10, "NormalisedH3K4LRorder.pdf")
draw(ht_list, merge_legend=T, heatmap_legend_side="bottom")
dev.off()




# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# Average over quantiles
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

q=10
quantile <- as.numeric(cut(1:length(H3K4_LRratio), q))


H3K4b_hP9C_profiles <- matrix(nrow=800, ncol=q)
H3K36_hP9C_profiles <- matrix(nrow=800, ncol=q)
hP9C_profiles <- matrix(nrow=800, ncol=q)
H3K4b_hP9C_profilesUN <- matrix(nrow=800, ncol=q)
H3K36_hP9C_profilesUN <- matrix(nrow=800, ncol=q)
hP9C_profilesUN <- matrix(nrow=800, ncol=q)
ZcwCVC_profiles <- matrix(nrow=800, ncol=q)
Zcw_hP9C_profilesUN <- matrix(nrow=800, ncol=q)
Input_hP9C_profiles <- matrix(nrow=800, ncol=q)

for(i in 1:q){
  H3K4b_hP9C_profiles[,i] <- colMeans(H3K4b_hP9C[ordering,][quantile==i,], na.rm = T) / colMeans(Input_hP9C[ordering,][quantile==i,], na.rm = T)
  H3K36_hP9C_profiles[,i] <- colMeans(H3K36_hP9C[ordering,][quantile==i,], na.rm = T) / colMeans(Input_hP9C[ordering,][quantile==i,], na.rm = T)
  hP9C_profiles[,i] <- colMeans(hP9C[ordering,][quantile==i,], na.rm = T) / colMeans(Input_hP9C[ordering,][quantile==i,], na.rm = T)
  ZcwCVC_profiles[,i] <- colMeans(Zcw_hP9C[ordering,][quantile==i,], na.rm = T) / colMeans(Zcw[ordering,][quantile==i,], na.rm = T)
  Zcw_hP9C_profilesUN[,i] <- colMeans(Zcw_hP9C[ordering,][quantile==i,], na.rm = T)
  H3K4b_hP9C_profilesUN[,i] <- colMeans(H3K4b_hP9C[ordering,][quantile==i,], na.rm = T)
  H3K36_hP9C_profilesUN[,i] <- colMeans(H3K36_hP9C[ordering,][quantile==i,], na.rm = T)
  hP9C_profilesUN[,i] <- colMeans(hP9C[ordering,][quantile==i,], na.rm = T)
  Input_hP9C_profiles[,i] <- colMeans(Input_hP9C[ordering,][quantile==i,], na.rm = T)
}

quantised_profiles <- rbind(data.table(H3K4b_hP9C_profiles, position=seq(-2000,1995,5), chip="Chip_H3K4me3_w/_hP9C vs In_hP9C"),
                            data.table(H3K36_hP9C_profiles, position=seq(-2000,1995,5), chip="Chip_H3K36me3_w/_hP9C vs In_hP9C"),
                            data.table(hP9C_profiles, position=seq(-2000,1995,5), chip="Chip_hP9C  vs In_hP9C"),
                            data.table(H3K4b_hP9C_profilesUN, position=seq(-2000,1995,5), chip="Chip_H3K4me3_w/_hP9C Unnormalised"),
                            data.table(H3K36_hP9C_profilesUN, position=seq(-2000,1995,5), chip="Chip_H3K36me3_w/_hP9C Unnormalised"),
                            data.table(hP9C_profilesUN, position=seq(-2000,1995,5), chip="Chip_hP9C Unnormalised"),
                            data.table(Input_hP9C_profiles, position=seq(-2000,1995,5), chip="Input_hP9C"),
                            data.table(ZcwCVC_profiles, position=seq(-2000,1995,5), chip="Zcw_hP9C vs Zcw"),
                            data.table(Zcw_hP9C_profilesUN, position=seq(-2000,1995,5), chip="Zcw_hP9C Unnormalised"))

pdf(width = 10, height = 10, file = "LR_Ratio_Profiles.pdf")
ggplot(melt(quantised_profiles, id.vars=c("position", "chip"), variable.name="H3K4LeftRightRatioDecile"),
       aes(position, value, colour=as.integer(H3K4LeftRightRatioDecile), group=H3K4LeftRightRatioDecile)) +
  geom_line() +
  scale_colour_gradient2(midpoint = q/2) +
  facet_wrap(~chip, scales = "free_y") +
  theme(legend.position = "bottom")

qplot(H3K4_LRratio, H3K36_LRratio, alpha=I(0.1)) +
  scale_y_log10() +
   scale_x_log10() +
   geom_abline(intercept = 0, slope = 1, col='red') +
   geom_smooth(method = "lm") +
   ggtitle("Input Normalised Left-Right Ratio\n of mean fragment depth H3K4 vs H3K36") +
   annotate("text", x = 0.3, y = 7, label = paste("Correlation:",round(cor(log(H3K4_LRratio), log(H3K36_LRratio)),3)))

qplot(H3K4_LRratioUN, H3K36_LRratioUN, alpha=I(0.1)) +
   scale_y_log10() +
   scale_x_log10() +
   geom_abline(intercept = 0, slope = 1, col='red') +
   geom_smooth(method = "lm") +
   ggtitle("UnNormalised Left-Right Ratio\n of mean fragment depth H3K4 vs H3K36") +
   annotate("text", x = 0.3, y = 7, label = paste("Correlation:",round(cor(log(H3K4_LRratioUN), log(H3K36_LRratioUN)),3)))

inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }
qplot(H3K4_LRratio, ZcwCVC_LRratio, alpha=I(0.1)) +
   scale_y_log10() +
   scale_x_log10() +
   geom_abline(intercept = 0, slope = 1, col='red') +
   geom_smooth(method = "lm") +
   ggtitle("Normalised Left-Right Ratio\n of mean fragment depth H3K4(/Input) vs Zcw_hP9C(/Zcw)") +
   annotate("text", x = 0.3, y = 7, label = paste("Correlation:",round(cor(inf2NA(log(H3K4_LRratio)),
                                                                           inf2NA(log(ZcwCVC_LRratio)), use="pairwise.complete.obs"),3)))

qplot(H3K4_LRratio, Zcw_hP9C_LRratioUN, alpha=I(0.1)) +
   scale_y_log10() +
   scale_x_log10() +
   geom_abline(intercept = 0, slope = 1, col='red') +
   geom_smooth(method = "lm") +
   ggtitle("UnNormalised Left-Right Ratio\n of mean fragment depth H3K4 vs Zcw_hP9C") +
   annotate("text", x = 0.3, y = 7, label = paste("Correlation:",round(cor(inf2NA(log(H3K4_LRratio)),
                                                                           inf2NA(log(Zcw_hP9C_LRratioUN)), use="pairwise.complete.obs"),3)))

dev.off()
