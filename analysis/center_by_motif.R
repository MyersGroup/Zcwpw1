

args=commandArgs(TRUE)

fastafile = args[1] #"motifs/HP9combo_peaks_NonPromoters_300w.fasta"
motiffile = args[2] #"motifs/Human_Motif_Results_Final_iter240.r"
inputbed = args[3] #"motifs/HP9N_peaks_NonPromoters_300w.bed"
outputbed = args[4] #"deeptools/beds/HP9N_peaks_NonPromoters_MotifCenteredStranded.bed"
plotfile = args[5] #"motifs/P9_motif_locations_at_P9.pdf"

# fastafile = "motifs/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL.fasta"
# motiffile = "motifs/Human_Motif_Results_Final_iter240.r"
# inputbed = "motifs/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_QCfiltered.bed"
# outputbed = "deeptools/beds/SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL_MotifCenteredStranded.bed"
# plotfile = "motifs/P9_motif_locations_SingleBasePeaks.SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140.p0.000001.sep250.ALL.pdf"

library(MotifFinder)
library(ggplot2)
library(data.table)

load(motiffile)

transform_published_motifs <- function(znew){
  compvec=c(0,1,0,0,0,0,1,1,0,0,1,1,0,0,0,1,1) #specify whether to take reverse complement of each motif for final output (0=no, 1=yes)
  keep=c(4,6,8,9,10,15,17,5,11,13,7,16,2,3,12,14,1) #specify the indices of motifs to keep (non-degenerate motifs)
  keep=c(4,6,8,9,10,15,17)


  prior=znew$prior
  starts=c(1,cumsum(znew$scorematdim)+1)
  starts=starts[1:length(znew$scorematdim)]
  ends=cumsum(znew$scorematdim)
  mot=matrix(nrow=0,ncol=4)
  sizes=vector(length=0)
  for(j in keep){
    scoremat = znew$scoremat[starts[j]:ends[j],]
    if(compvec[j]==1){
      compmat=scoremat[,c(4:1)]
      compmat=compmat[nrow(compmat):1,]
      scoremat=compmat
    }
    mot=rbind(mot,scoremat)
    sizes=c(sizes,nrow(scoremat))
  }
  return(list(mot=mot, sizes=sizes, prior=prior))
}

znew <- transform_published_motifs(znew)

# Plot motifs to check
# tmp <- t(exp(znew[[1]]))
# rownames(tmp) <- c("A", "C", "G", "T")
#
# pdf(width=10, height=5, file="motifs/Altemose_Prdm9.pdf")
# pos <- cumsum(c(1,znew[[2]]))
# for(i in 1:7){
#   print(ggseqlogo::ggseqlogo(tmp[,pos[i]:(pos[i+1]-1)]))
# }
# dev.off()

HP9N_peaks_NonPromoters_300w <- load_sequences(fastafile, sm.rm = F)

HP9N_peaks_NonPromoters_300w_positions <- getmotifs(znew$mot, znew$sizes,
                                                    HP9N_peaks_NonPromoters_300w, maxwidth = max(nchar(HP9N_peaks_NonPromoters_300w)),
                                                        alpha = 0.2, maxits = 10, updatemot = 0, ourprior = znew$prior, seed=42, stranded_prior = T)

saveRDS(HP9N_peaks_NonPromoters_300w_positions, paste0(inputbed,"_inferredMotifPositions.rds"))
#HP9N_peaks_NonPromoters_300w_positions <- readRDS(paste0(inputbed,"_inferredMotifPositions.rds"))

export_locations <- function(HP9N_peaks_NonPromoters_300w_positions){
  tmp <- HP9N_peaks_NonPromoters_300w_positions$dt[!is.na(whichpos)]
  tmp$sequence <- as.integer(tmp$sequence)
  setorder(tmp, sequence)

  motifcenters <- data.table(whichmotif=1:7,
                             motifcenterZF7C=c(19, 21, 21, 26, 11, 11, 19),
                             motifcenterZF4T=c(11, 12, 12, 21, 6, 6, 14), # NB motifcenterZF4T makes no disernable difference
                             motif_len=HP9N_peaks_NonPromoters_300w_positions$scorematdim)

  # head(extract_matches(HP9N_peaks_NonPromoters_300w_positions))

  tmp <- tmp[motifcenters, on="whichmotif"]
  tmp[whichstrand==0, whichstrand := -1]

  tmp[whichstrand==1, whichpos := whichpos + motifcenterZF7C]
  tmp[whichstrand==-1,whichpos := whichpos + motif_len - motifcenterZF7C]

  tmp[,whichpos := whichpos-150] # assuming slop of 150 in rule extractFASTA

  tmp[whichstrand==-1, strand := "-"]
  tmp[whichstrand==1, strand := "+"]

  return(tmp)
}

tmp <- export_locations(HP9N_peaks_NonPromoters_300w_positions)

HP9N_peaks_NonPromoters <- fread(inputbed)
HP9N_peaks_NonPromoters$sequence <- 1:nrow(HP9N_peaks_NonPromoters)

# Diagnosis plots
pdf(width=10, height=5, file=plotfile)
plot_motif_location(HP9N_peaks_NonPromoters_300w_positions, top_n=1000)

plot(HP9N_peaks_NonPromoters_300w_positions$prior)

# check convergence
ggplot(melt(data.table(cbind(HP9N_peaks_NonPromoters_300w_positions$alphas, seq=1:10)), id="seq"), aes(seq, value, colour=variable)) + geom_line()

#Regprob_by_enrichment
plot(HP9N_peaks_NonPromoters[tmp, on="sequence"][order(sequence)]$regprob, pch=".")
dev.off()

# v1 = chr
# v9 = enrichment
# v10 = likelihood

fwrite(HP9N_peaks_NonPromoters[tmp, on="sequence"][,.(V1, V2+whichpos, V3+whichpos, sequence, V9, strand, whichmotif, V10, V6, V7, V8)][order(sequence)], file=outputbed, col.names = FALSE, sep = "\t")



