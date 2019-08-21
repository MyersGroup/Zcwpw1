
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
