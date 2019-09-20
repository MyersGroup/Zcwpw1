library(data.table)
library(ggplot2)
library(purrr)
library(dplyr)

readdmc1 <- function(filename){
  sample <- fread(filename)

  if(grepl("5prime",filename)){
    sample$Strand <- "5 prime"
  }else{
    sample$Strand <- "3 prime"
  }

  if(grepl("atKO",filename)){
    sample$Genotype <- "B6 Prdm9 KO DMC1 regions"
  }else{
    sample$Genotype <- "B6 Wild Type DMC1 regions"
  }

  if(grepl("Brick2012",filename)){
    sample$Mouse <- "WT"
  }else{
    sample$Mouse <- "Zcwpw1-/-"
  }

  return(sample)
}

tmp <-
  list.files(pattern = "*prime_at*") %>%
  map_df(~readdmc1(.))

names(tmp)[1:2] <- c("Position","MeanCoverage")
tmp[, Genotype := factor(Genotype, levels=c("B6 Wild Type DMC1 regions","B6 Prdm9 KO DMC1 regions"))]

tmp[Position<(-3000), background := mean(MeanCoverage), by=c("Strand","Mouse")]
tmp[, background := mean(background, na.rm = T), by=c("Strand","Mouse")]

tmp[, MeanCoverage := MeanCoverage - background, by=c("Strand","Mouse")]

tmp[, MeanCoverage := MeanCoverage / sum(MeanCoverage), by=c("Strand","Mouse")]

pdf("DMC1_SSDS.pdf", height = 5, width = 8)
ggplot(tmp[Position>(-2000)], aes(Position, MeanCoverage, colour=Mouse, alpha=Strand)) + #group=Strand,
  geom_line(size=0.35) +
  facet_wrap(~Genotype, nrow=1) +
  scale_color_brewer(palette = "Set1") +
  xlab("Position from center (bp)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.border = element_rect(color = "grey40", fill = NA, size = 0.5),
        strip.background = element_rect(color = "grey40", size = 0.5)) +
  scale_alpha_manual(values=c(1,0.5)) +
  ylab("Normalised Mean Coverage")
dev.off()
