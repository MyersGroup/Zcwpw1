library(data.table)
library(ggplot2)

args = commandArgs(TRUE)

outfile <- args[1]
regionsName <- args[2]

plotwidth <- as.numeric(args[3])
plotheight <- as.numeric(args[4])

args <- args[5:length(args)]

# test values
# args <- c("bwprofilesNorm/SRA_Altemose2015_SRR5627149_VS_SRA_Altemose2015_SRR5627143_AT_SingleBasePeaks.SRA_Altemose2015_SRR5627146_AND_SRA_Altemose2015_SRR5627147_vs_SRA_Altemose2015_SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded.profile",
#           "bwprofilesNorm/SRA_Altemose2015_SRR5627146_AND_SRA_Altemose2015_SRR5627147_VS_SRA_Altemose2015_SRR5627143_AT_SingleBasePeaks.SRA_Altemose2015_SRR5627146_AND_SRA_Altemose2015_SRR5627147_vs_SRA_Altemose2015_SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded.profile",
#           "ChipH3K36_hP9HA vs In_hP9HA",
#           "ChipHA_hP9HA+V5 vs In_hP9HA")

files <- args[1:(length(args)/2)]
names <- args[(length(args)/2+1):length(args)]

enrichment <- data.table(Position=integer(),
                      Q1=numeric(),
                      Q2=numeric(),
                      Q3=numeric(),
                      Q4=numeric(),
                      Random=numeric(),
                      samplePair=character())

for(i in seq_along(files)){
  profile <- fread(files[i])
  setnames(profile, c("Position","Q1","Q2","Q3","Q4","Random"))
  profile[,samplePair := names[i]] #gsub("bwprofilesNorm/","",strsplit(file,"_AT_")[[1]][1])
  enrichment <- rbind(enrichment, profile)
}


enrichment[, samplePair := gsub(" vs","\n vs",samplePair)]

pdf(outfile, width = plotwidth, height = plotheight)
ggplot(melt(enrichment, id.vars=c("Position","samplePair"), variable.name="Region", value.name="Enrichment"), aes(Position, Enrichment, colour=Region)) +
  geom_hline(yintercept = 1) +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  scale_y_log10() +
  facet_wrap(~samplePair, scales='free_y') +
  expand_limits(y = c(1, 2)) +
  ggtitle(regionsName) #gsub(".pdf","",strsplit(outfile,"_AT_")[[1]][2])
dev.off()
