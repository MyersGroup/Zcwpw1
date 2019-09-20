

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

# generate profile data
snakemake --cores 15 --snakefile pipelines/Plot_Heatmap.py -npr

library(data.table)
library(ComplexHeatmap)


normalise_matrix <- function(sample, input, pseudocount=1, range=c(0.025, 1-0.025)){

  norm <- (sample+pseudocount)/(input+pseudocount)

  # 0/0 = 0
  norm[is.na(norm)] <- 0

  # x/0 = Inf, pass
  norm[!is.finite(norm)] <- NA

  # Clip extreme values
  norm[norm>quantile(norm, na.rm=T, range[2])] <- quantile(norm, na.rm=T, range[2])
  norm[norm<quantile(norm, na.rm=T, range[1])] <- quantile(norm, na.rm=T, range[1])

  return(norm)
}




makeColumn <- function(sample, input, title, ordering){

  # create line plot for top of heatmap
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


read_normmatrix <- function(sample, locations, const){
  # read in matrix, normalised by total read count
  # const to avoid numeric instability
  return(
    as.matrix(fread(paste0("bwmatrices/",sample,"_AT_",locations,".bwm"))) /
      (as.numeric(fread(paste0("FragPos/Fragment_Position_",sample,".total"))$V1)/const)
  )
}


plotHeatmap <- function(locations, samples, controls, names, title){

  const <- 1e10

  # first samples are used to create ordering only
  tmp <- normalise_matrix(read_normmatrix(samples[1], locations, const),
                          read_normmatrix(controls[1], locations, const))
  ordering <- order(-apply(tmp[,(ncol(tmp)/2-100):(ncol(tmp)/2+100)], 1, function(x) mean(x, na.rm =T)))

  # create list of heatmap column objects
  for(i in seq_along(samples)[-1]){
    if(i==2){ # create list
      ht_list <- makeColumn(read_normmatrix(samples[i], locations, const),
                            read_normmatrix(controls[i], locations, const),
                            names[i],
                            ordering)
    }else{ # append to list
      ht_list <- ht_list + makeColumn(read_normmatrix(samples[i], locations, const),
                                      read_normmatrix(controls[i], locations, const),
                                      names[i],
                                      ordering)
    }
  }

  pdf(width = 16, height = 10, paste0("bwplots/NormalisedhP9C_vInOrder_",title,"-locations.pdf"))
  draw(ht_list, merge_legend=T, heatmap_legend_side="bottom")
  dev.off()

}

# first sample is used to set up the row ordering only

samples <- c("WTCHG_538916_223180",
           #  "WTCHG_538916_221156",
          #   "WTCHG_538916_223180",
             "WTCHG_538916_223180",
          "WTCHG_538916_224192",
             "NA15-SRR5627152_AND_NA15-SRR5627153",
             "NA15-SRR5627149",
             "NA15-SRR5627146_AND_NA15-SRR5627147",
          "NA15-SRR5627145_AND_NA15-SRR5627144",
             "NA15-SRR5627150_AND_NA15-SRR5627151",
             "NA15-SRR5627148")

controls <- c("WTCHG_538916_221156",
             # "WTCHG_538916_217108_AND_WTCHG_538916_220144",
            #  "WTCHG_538916_217108_AND_WTCHG_538916_220144",
              "WTCHG_538916_221156",
            "WTCHG_538916_221156",
              "NA15-SRR5627143",
              "NA15-SRR5627143",
              "NA15-SRR5627143",
            "NA15-SRR5627143",
              "NA15-SRR5627142",
              "NA15-SRR5627142")

names <- c("ZHA_hP9V5 vs ZHA",
           #"ZHA vs In",
           #"ZHA_hP9V5 vs In",
           "Zcwpw1 with hPrdm9 vs w\\o", #"ZHA_hP9V5 vs ZHA",
           "Zcwpw1 with cPrdm9 vs w\\o", #"ZHA_hP9V5 vs ZHA",
           "H3K4me3 with hPrdm9", #"H3K4_hP9HA vs InhP9HA",
           "H3K36me3 with hPrdm9", #"H3K36_hP9HA vs InhP9HA",
           "hPrdm9", #"hP9HA vs InhP9HA",
           "cPrdm9", #"hP9HA vs InhP9HA",
           "Untransfected H3K4",
           "Untransfected H3K36")

plotHeatmap("SingleBasePeaks.WTCHG_538916_223180_vs_WTCHG_538916_221156.p0.000001.sep250.ALL_Q00", samples, controls, names, "hZcwCTvsST")
plotHeatmap("SingleBasePeaks.WTCHG_538916_221156_vs_WTCHG_538916_217108.p0.000001.sep250.ALL_Q00", samples, controls, names, "hZcwSTvsIn")
plotHeatmap("SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded_Q00", samples, controls, names, "HP9C")
plotHeatmap("SingleBasePeaks.NA15-SRR5627138_AND_NA15-SRR5627139_vs_NA15-SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded_Q00", samples, controls, names, "HP9N")

plotHeatmap("SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL_Q00", samples, controls, names, "CHP9C")

