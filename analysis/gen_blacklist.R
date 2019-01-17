
args <- commandArgs(TRUE)

str(args)

input <- args[1]
output <- args[2]

library(data.table)
library(ggplot2)
cov_ic <- fread(paste0("zcat < ",input))
names(cov_ic) <- make.names(names(cov_ic))

names(cov_ic)[4] <- "sample"

cov_ic <- cov_ic[X.chrom!="chrM"] # only one data point

pdf(file = paste0(output,"_coverage.pdf"))

ggplot(cov_ic, aes(start, sample)) +
        geom_line() +
        scale_y_sqrt() +
		geom_hline(yintercept=quantile(cov_ic$sample, 0.9999), colour='red') +
        facet_wrap(~X.chrom, scales="free_x")

ggplot(cov_ic, aes(start, sample)) +
	geom_line() +
	coord_cartesian(ylim=c(0,max(cov_ic$sample))) +
	geom_hline(yintercept=quantile(cov_ic$sample, 0.9999), colour='red') +
	facet_wrap(~X.chrom, scales="free_x")

ggplot(cov_ic, aes(start, sample)) +
	geom_line() +
	coord_cartesian(ylim=c(0,quantile(cov_ic$sample, 0.99999))) +
	geom_hline(yintercept=quantile(cov_ic$sample, 0.9999), colour='red') +
	facet_wrap(~X.chrom, scales="free_x")

ggplot(cov_ic, aes(start, sample)) +
	geom_line() +
	coord_cartesian(ylim=c(0,quantile(cov_ic$sample, 0.9999))) +
	facet_wrap(~X.chrom, scales="free_x")

ggplot(cov_ic[X.chrom %in% paste0("chr",1) & start!=0], aes(start, sample)) +
        geom_line() +
		geom_hline(yintercept=quantile(cov_ic$sample, 0.9999), colour='red') +
        coord_cartesian(ylim=c(0,max(cov_ic[X.chrom %in% paste0("chr",1) & start!=0]$sample)))
dev.off()

print(paste("Using blacklist threshold value of", quantile(cov_ic$sample, 0.9999)))

write.table(cov_ic[sample>quantile(cov_ic$sample, 0.9999)][,.(X.chrom, start, end)],
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t",
            file = output)
