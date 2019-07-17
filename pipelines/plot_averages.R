
# used for byfragmentlength analysis

library(data.table)
library(ggplot2)
library(parallel)

args = commandArgs(TRUE)

outfile <- args[1]

args <- args[2:length(args)]

# Args should be a list of files, first x being samples, second x being matched input samples

# args = c("byfraglength/WTCHG_538916_223180_00.profile","byfraglength/WTCHG_538916_223180_01.profile",
#         "byfraglength/WTCHG_538916_223180_02.profile","byfraglength/WTCHG_538916_223180_03.profile", "byfraglength/WTCHG_538916_223180_04.profile",
#         "byfraglength/WTCHG_538916_221156_00.profile","byfraglength/WTCHG_538916_221156_01.profile",
#         "byfraglength/WTCHG_538916_221156_02.profile","byfraglength/WTCHG_538916_221156_03.profile", "byfraglength/WTCHG_538916_221156_04.profile")


args <- gsub("\\.[a-z]*","",args)

if(length(args)>24){stop("Too Many Files, Maximum is 12 sample-input pairs (24 total)")}


read_freq <- function(file){
  data.table(fread(file, integer64 = "double"), sample=gsub("\\.[a-z]*","",file))
}

read_and_normalise <- function(args){
  numCores <- detectCores()

  freqs <- mclapply(paste0(args,".profile"), read_freq, mc.cores = max(c(12,numCores)))

  freqs <- rbindlist(freqs)

  # total number of reads to normalise
  totals <- mclapply(paste0(args,".total"), read_freq, mc.cores = max(c(12,numCores)))
  totals <- rbindlist(totals)

  setkey(freqs, sample)
  setkey(totals, sample)

  tmp <- merge(freqs, totals)

  setnames(tmp,"V1.x","Position")
  tmp[,Mean_CPM := V2 / (V1.y/1e6)]

  # assuming file has been split into quantiles which are the last two characters of the file
  tmp[, Quantile := stringr::str_sub(sample, start= -2)]

  setkey(tmp, Position, Quantile)

  return(tmp)
}

sample <- read_and_normalise(args[1:(length(args)/2)])
input <- read_and_normalise(args[(length(args)/2 + 1):length(args)])

# tmp <- merge(freqs, totals)
#
# setnames(tmp,"V1.x","Position")
# tmp[,Mean_CPM := V2 / (V1.y/1e6)]
#
# tmp2 <- merge(freqs, totals)
# setnames(tmp2,"V1.x","Position")
# tmp2[,Mean_CPM := V2 / (V1.y/1e6)]
#
# tmp[, Quantile := stringr::str_sub(sample, start= -2)]
# tmp2[, Quantile := stringr::str_sub(sample, start= -2)]
#
# setkey(tmp, Position, Quantile)
# setkey(tmp2, Position, Quantile)


pdf(outfile)
ggplot(merge(sample, input), aes(Position, Mean_CPM.x/Mean_CPM.y, colour=Quantile)) +
  geom_line() +
  xlab("Position from Motif Centre") +
  ylab("Average Fragment Depth per total fragment coverage") +
  ggtitle("Zcw_hP9C Co-transfection normalised by Zcw alone (ChipVsChip)") +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom")
dev.off()

