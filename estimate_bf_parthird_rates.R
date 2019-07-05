
library(reshape2)

######### final cosine similarity plots & Neff (obs/sim parallel in diploid)
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

# parsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", as.is = T)
INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_1plus1_only_writefrac/"
i1files <- list.files(path = INDIR, pattern = "_infsites_permut_effgenfrac0.1.txt", full.names = T, recursive = T)

i1samples <- sub(pattern = "_infsites_permut_effgenfrac0.1.txt", replacement = "", x = basename(i1files))
i1df <- as.data.frame(do.call(rbind, lapply(X = i1files, FUN = function(x) read.delim(file = x, as.is = T)$counts_tot)))
# i1df$sampleid <- i1samples

reftable <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isavioltypes_reffile.txt", as.is = T)

i1df <- as.data.frame(t(apply(X = i1df, MARGIN = 1, FUN = function(x, f) by(data = x, INDICES = f, FUN = sum), f = reftable$type)))
i1df$sampleid <- i1samples

# i1df <- reshape(data = i1df, direction = "wide", idvar = "sampleid", timevar = "type")

i1df$totviol <- rowSums(i1df[, 1:4])
i1df$backrate <- i1df$back/i1df$totviol
i1df$forwardrate <- i1df$forward/i1df$totviol
i1df$parrate <- i1df$parallel/i1df$totviol
i1df$thirdrate <- i1df$third_allele/i1df$totviol

i1df$histology_abbreviation <- sumtable_whole[match(x = i1df$sampleid, table = sumtable_whole$samplename), "histology_abbreviation"]
write.table(x = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/bf_pt_rate_estimates.txt", quote = F, row.names = F, sep = "\t")

parsamplesdfm <- merge(x = parsamplesdf, y = i1df, by = "sampleid", all.x = T, all.y = F)

library(ggplot2)

plotdf <- parsamplesdfm[order(parsamplesdfm$med_diploid, decreasing = T), ]
plotdf$sampleid <- factor(x = plotdf$sampleid, levels = plotdf$sampleid)
plotdf <- plotdf[which(plotdf$tot_diploid/plotdf$tot_testable > .10 & plotdf$freq_med.parallel >= 10), ]
plotdf$neff <- plotdf$med_diploid/(plotdf$freq_med.parallel/10)