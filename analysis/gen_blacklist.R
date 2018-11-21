
args <- commandArgs(TRUE)

print(str(args))

bedloc <- args[1]
id <- args[2]
output_dir <- args[3]

library(data.table)
library(ggplot2)
cov_ic <- fread(paste0("zcat < ",bedloc,id,"_indexcov-indexcov.bed.gz"))
names(cov_ic) <- make.names(names(cov_ic))

names(cov_ic)[4] <- "sample"

pdf(file = paste0(output_dir,id,"_coverage.pdf"))

#[X.chrom %in% paste0("chr",1:22) & start!=0]
ggplot(cov_ic, aes(start, sample)) +
        geom_line() +
        scale_y_sqrt() +
        facet_wrap(~X.chrom, scales="free_x")

ggplot(cov_ic, aes(start, sample)) +
        geom_line() +
        coord_cartesian(ylim=c(0,1000)) +
        facet_wrap(~X.chrom, scales="free_x")

ggplot(cov_ic[X.chrom %in% paste0("chr",10) & start!=0], aes(start, sample)) +
        geom_line() +
        coord_cartesian(ylim=c(0,1000))
dev.off()

write.table(cov_ic[sample>10][,.(X.chrom, start, end)],
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t",
            file = paste0(id,"_blacklist.bed"))
