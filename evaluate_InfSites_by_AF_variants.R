#### evaluate parallel calls

library(ggplot2)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

allsumfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", pattern = "_InfSites_VAFpipeline_summarystats.txt", full.names = T, recursive = T)
allsums <- lapply(X = allsumfiles, FUN = read.table, row.names = 1, as.is = T)
allsums <- data.frame(t(do.call(cbind, allsums)))
allsums$sampleid <- gsub(pattern = "_InfSites_VAFpipeline_summarystats.txt", replacement = "", x = basename(allsumfiles))

allpcfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
allpc <- lapply(X = allpcfiles, FUN = read.table, col.names = c("pseudocount", "snv_slope"), as.is = T)
allpc <- data.frame(do.call(rbind, allpc))
allpc$sampleid <- gsub(pattern = "_pseudocount_calibrated.txt", replacement = "", x = basename(allpcfiles))

allsums <- merge(x = allsums, y = allpc, by = "sampleid")

sumtab <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)
allsums$histology_abbreviation <- sumtab[match(x = allsums$sampleid, table = sumtab$samplename), "histology_abbreviation"]

#this is THE set to work with:
View(allsums[which(allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1), ])

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(allsums$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(allsums$histology_abbreviation)

p1 <- ggplot(data = allsums[allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1, ], mapping = aes(x = med, y = nparallel, colour = histology_abbreviation)) + geom_point()
p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks() + geom_abline()
p1

p1 <- ggplot(data = allsums[allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1, ], mapping = aes(x = tot_hetero, y = nparallel, colour = npar_phased>0)) + geom_point()
p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks()
p1

p1 <- ggplot(data = allsums, mapping = aes(x = tot_hetero, y = nparallel, colour = log2(nparallel/tot_hetero)>-8)) + geom_point()
p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks()
p1

p1 <- ggplot(data = allsums, mapping = aes(x = log2(nparallel/tot_hetero), fill = npar_phased>0)) + geom_histogram()
p1


write.table(x = allsums, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary.txt", quote = F, col.names = T, row.names = F)
write.table(x = allsums[which(allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1), ], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", quote = F, col.names = T, row.names = F)
