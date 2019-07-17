library(data.table)
library(ggplot2)
library(parallel)

args = commandArgs(TRUE)

outfile <- args[1]

args <- args[2:length(args)]

# args <- c("bwprofiles/SRA_Altemose2015_SRR5627146_AND_SRA_Altemose2015_SRR5627147_AT_SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.profile",
#           "14504698649",
#           "bwprofiles/SRA_Altemose2015_SRR5627143_AT_SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.profile",
#           "9095954710")

sample <- fread(args[1])
control <- fread(args[3])

sample <- as.matrix(sample)
control <- as.matrix(control)

enrichment <- (sample[,2:ncol(sample)] / (as.numeric(args[2])/1e6) ) / (control[,2:ncol(control)] / (as.numeric(args[4])/1e6) )

enrichment <- data.table(Position=sample[,1], enrichment)

setnames(enrichment, c("Position","Q1","Q2","Q3","Q4","Random"))

pdf(outfile)
ggplot(melt(enrichment, id.vars="Position", variable.name="Quantile", value.name="Enrichment"), aes(Position, Enrichment, colour=Quantile)) + geom_line()
dev.off()
