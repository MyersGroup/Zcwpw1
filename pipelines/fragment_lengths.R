
# Rscript fragment_lengths.R fraglen.pdf MAPeakCaller/Fragment_Position_SRA_Altemose2015_SRR5627152.sorted.bed MAPeakCaller/Fragment_Position_SRA_Altemose2015_SRR5627143.sorted.bed

library(data.table)
library(ggplot2)
library(parallel)

args = commandArgs(TRUE)

outfile <- args[1]

args <- args[2:length(args)]

if(length(args)>12){stop("Too Many Files, Maximum is 12")}

#args = c("MAPeakCaller/Fragment_Position_SRA_Altemose2015_SRR5627152.sorted.bed","MAPeakCaller/Fragment_Position_SRA_Altemose2015_SRR5627143.sorted.bed")

read_freq <- function(file){
  data.table(fread(paste0("awk '{print $3-$2; }' ",file))[,.N,by=V1][order(-V1)], sample=file)
}

numCores <- detectCores()

freqs <- mclapply(args, read_freq, mc.cores = max(c(12,numCores)))

frequencies <- rbindlist(freqs)

frequencies[ , percentage := prop.table(N) , by = "sample"]

pdf(outfile)
ggplot(frequencies, aes(V1, percentage, colour=sample)) +
  geom_line() +
  xlab("Fragment length") +
  scale_colour_brewer(palette = "Paired") +
  theme(legend.position = "bottom") +
  guides(colour=guide_legend(ncol=1)) +
  xlim(0,800)
dev.off()
