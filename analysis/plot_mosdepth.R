
# e.g. to run this script
# Rscript plot_mosdepth.R "594404" #"278163|276139|217108|220144"

args <- commandArgs(TRUE)

print(str(args))

project <- args[1]

library(data.table)
library(ggplot2)
library(gghighlight)


files=grep(project, list.files("qc/mosdepth/", pattern=".txt"), value = T)
name=project

        coverage <- data.table()
        for (i in files){
                tmp <- fread(paste0("qc/mosdepth/", i))
                tmp$sample <- i
                coverage <- rbind(coverage, tmp)
        }
        
        coverage[,Stage:="De-Duplicated"]
        coverage[grepl("raw",sample),Stage:="Raw"]
        coverage[, sample:=gsub("raw_|.mosdepth.global.dist.txt","",perl = T, sample)]
		coverage[, sample:=gsub(name,"",perl = T, sample)]
        
        p <- ggplot(coverage[V1=="total"], aes(V2, V3, colour=Stage)) +
                geom_line() +
                scale_x_continuous(breaks=0:10, limits=c(0,10), minor_breaks = NULL)+
                scale_y_continuous(labels = scales::percent) +
                facet_wrap(~sample) +
                labs(x="Coverage", y="% of bases with at least >= this level of coverage") +
                theme(legend.position = "bottom", ) + 
                scale_color_brewer(palette="Set1") +
                guides(colour=guide_legend(nrow=2,byrow=TRUE))
        
pdf(height=10, width=15, file=paste0("qc/mosdepth/mosdepth_",name,".pdf"))

p

decoys <- c("chrUn_KN707896v1_decoy","chrUn_KI270438v1","chr14_GL000225v1_random","chrUn_KI270442v1","chr22_KI270736v1_random","chrUn_GL000216v2")

coverage[,chr:=as.integer(gsub("chr","",V1))]


# coverage by chrom over multitple %
ggplot(coverage[Stage=="De-Duplicated" & V1=="chrM"], aes(V2+1, V3, colour=sample)) +
  geom_line() +
  scale_x_log10() +
  gghighlight(use_direct_label=FALSE) +
  ggtitle("ChrM") +
  facet_wrap(~sample)+
  labs(x="Coverage +1", y="% of bases with at least >= this level of coverage")

  ggplot(coverage[Stage=="De-Duplicated" & V1=="chr1" &  V3>0.005], aes(V2+1, V3, colour=sample)) +
    geom_line() +
    scale_x_log10() +
    gghighlight(use_direct_label=FALSE) +
    ggtitle("Chr1") +
    facet_wrap(~sample)+
    labs(x="Coverage +1", y="% of bases with at least >= this level of coverage")



#ggplot(coverage[(nchar(V1)<=5 | V1 %in% decoys) & Stage=="De-Duplicated" &  V3>0.005], aes(V2+1, V3, colour=sample, linetype=sample)) +
#  geom_line() +
#  facet_wrap(~V1, scales='free_x') +
#  scale_x_log10() +
#  labs(x="Coverage +1", y="% of bases with at least >= this level of coverage")

  # coverage by chrom
#  ggplot(coverage[nchar(V1)<=5 & Stage=="De-Duplicated" & V2==3 & !is.na(chr)], aes(chr, V3)) + geom_col() + facet_wrap(~sample)


#  for(chrom in c(unique(coverage$V1)[1:25],decoys)){
#print(	  ggplot(coverage[Stage=="De-Duplicated" & V1==chrom], aes(V2+1, V3, colour=sample)) +
#	    geom_line() +
#	    scale_x_log10() +
#	    gghighlight(use_direct_label=FALSE) +
#	    ggtitle(chrom) +
#	    facet_wrap(~sample)+
#	    labs(x="Coverage +1", y="% of bases with at least >= this level of coverage")
#		)
 # 	
  #}

dev.off()