
# e.g. to run this script
# Rscript plot_mosdepth.R "594404" "278163|276139|217108|220144"

args <- commandArgs(TRUE)

print(str(args))

project <- args[1]
input_regex <- args[2]

library(data.table)
library(ggplot2)

plot_depth <- function(files, name="", input=input_regex){
        coverage <- data.table()
        for (i in files){
                tmp <- fread(paste0("mosdepth/", i))
                tmp$sample <- i
                coverage <- rbind(coverage, tmp)
        }
        
        coverage[,raw:="De-Duplicated"]
        coverage[grepl("raw",sample),raw:="Raw"]
        coverage[, sample:=gsub("raw_|.mosdepth.global.dist.txt","",perl = T, sample)]
        coverage[grepl(input,sample),sample:=paste("Input",sample)]
        coverage[!grepl(input,sample),sample:=paste("Chip",sample)] #"277151|279175"
        str(coverage)
        
        p <- ggplot(coverage[V1=="total"], aes(V2, V3, colour=sample)) +
                geom_line() +
                scale_x_continuous(breaks=0:10, limits=c(0,10), minor_breaks = NULL)+
                scale_y_continuous(labels = scales::percent) +
                facet_wrap(~raw, nrow=2) +
                labs(x="Coverage", y="% of bases with at least >= this level of coverage") +
                theme(legend.position = "bottom", ) + 
                scale_color_brewer(palette="Set1") +
                guides(colour=guide_legend(nrow=2,byrow=TRUE))
        
        ggsave(paste0("mosdepth/mosdepth_",name,".pdf"), p, width = 8, height=8)
        print(p)
        
        return(coverage[V1=="total"][order(V2)])
}

plot_depth(grep(project, list.files("mosdepth/"), value = T), name=project)

