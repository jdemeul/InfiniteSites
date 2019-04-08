### evaluating fraction of SNVs which may be picked up based on power

library(parallel)

parfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/", pattern = "_baflogr_gc_results.txt", recursive = T, full.names = T)
parmuts <- mclapply(X = parfiles, FUN = read.delim, as.is = T, mc.cores = 16, mc.preschedule = T)
