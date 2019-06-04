# get samples for which this has completed

INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_het_only_writefrac/"
p1files <- list.files(path = INDIR, pattern = "_infsites_totals_effgenfrac0.1.txt", full.names = T, recursive = T)
i1files <- list.files(path = INDIR, pattern = "_infsites_totals_effgenfrac1.txt", full.names = T, recursive = T)

p1samples <- sub(pattern = "_infsites_totals_effgenfrac0.1.txt", replacement = "", x = basename(p1files))
i1samples <- sub(pattern = "_infsites_totals_effgenfrac1.txt", replacement = "", x = basename(i1files))

commonsamples <- intersect(p1samples, i1samples)

comparisons <- lapply(X = commonsamples, FUN = function(sampleid) {
  p1file <- grep(value = T, pattern = sampleid, x = p1files)
  i1file <- grep(value = T, pattern = sampleid, x = i1files)

  p1 <- read.delim(p1file, as.is = T)[3, 2:4]
  i1 <- read.delim(i1file, as.is = T)[3, 2:4]
  
  return(c(as.numeric(p1), as.numeric(i1)))
})

comparisons <- as.data.frame(do.call(rbind, comparisons))
colnames(comparisons) <- c("freq_low_p1", "freq_med_p1", "freq_hi_p1", "freq_low_1", "freq_med_1", "freq_hi_1")
comparisons$sampleid <- commonsamples

library(ggplot2)
library(reshape2)

plotdf <- melt(comparisons, id.vars = "sampleid", measure.vars = c("freq_med_p1", "freq_med_1"))
plotdf$f <- ifelse(plotdf$variable == "freq_med_p1", 0.1, 1)
pl1 <- ggplot(data = plotdf, mapping = aes(x = f, y = value)) + geom_point() + geom_path(mapping = aes(group = sampleid)) + scale_x_log10() + scale_y_log10() + annotation_logticks()
pl1

comparisons$slope <- (log10(comparisons$freq_med_p1)-log10(comparisons$freq_med_1))/(-1-0)
