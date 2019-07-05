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
consrhopsi <- read.delim(file = "/srv/shared/vanloo/ICGC_consensus_copynumber/consensus.20170217.purity.ploidy.txt.gz", as.is = T)
allsums <- cbind(allsums, consrhopsi[match(x = allsums$sampleid, table = consrhopsi$samplename), -1])
allsums$is_preferred <- allsums$sampleid %in% sumtab[sumtab$is_preferred, "samplename"]

# allsums$hitratio <- log10(allsums$nparallel/allsums$tot_testable)
# allsums$hitratio2 <- log10(allsums$med/allsums$tot_testable)
allsums$hitratio3 <- log10(allsums$tot_testable^2/allsums$nparallel)
allsums$hitratio4 <- log10(allsums$tot_testable^2/allsums$med)
# allsums$hitratio5 <- log10(allsums$tot_diploid^2/allsums$med_diploid)
allsums$pratio <- pnorm(q = allsums$hitratio3, mean = mean(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), sd = sd(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), lower.tail = TRUE, log.p = FALSE)

# p1 <- ggplot(data = allsums[which(allsums$prec > 0 & allsums$rec > 0), ]) + geom_point(mapping = aes(x = hitratio2, y = purity_conf_mad))
# p1

# p1 <- ggplot(data = allsums) + geom_point(mapping = aes(x = hitratio, y = purity_conf_mad, colour = hitratio > -2.352241))
# p1

# ratecutoff <- qnorm(p = .01, mean = mean(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), sd = sd(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), lower.tail = TRUE, log.p = FALSE)
p1 <- ggplot(data = allsums) + geom_point(mapping = aes(y = log10(nparallel), x = log10(tot_testable^2), colour = pratio > .01))
p1

p1 <- ggplot(data = allsums) + geom_histogram(mapping = aes(x = hitratio3, fill = !is.na(allsums$prec) & allsums$prec > 0.5 & allsums$rec > 0))
p1

p1 <- ggplot() + geom_histogram(data = allsums[allsums$tot_phaseable >= 1000, ], mapping = aes(x = hitratio4), fill = "blue", alpha = .5) +
  geom_histogram(data = allsums, mapping = aes(x = hitratio3), fill = "red", alpha = .5)
p1

p1 <- ggplot() + geom_histogram(data = allsums[allsums$tot_phaseable >= 1000, ], mapping = aes(x = hitratio4), fill = "blue", alpha = .5) +
  geom_histogram(data = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], mapping = aes(x = hitratio3), fill = "red", alpha = .5)
p1


# p1 <- ggplot() + geom_point(data = allsums[allsums$tot_phaseable > 1000, ], mapping = aes(y = log10(med), x = log10(tot_testable^2)), colour = "black") +
#   geom_point(data = allsums, mapping = aes(y = log10(nparallel), x = log10(tot_testable^2), colour = hitratio > -2.3))
# p1



#this is THE set to work with:
View(allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ]) #ALT
# View(allsums[which(allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1), ])

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(allsums$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(allsums$histology_abbreviation)

p1 <- ggplot(data = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], mapping = aes(x = med, y = nparallel, colour = histology_abbreviation)) + geom_point()
p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks() + geom_abline()
p1

p1 <- ggplot(data = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], mapping = aes(x = tot_hetero, y = nparallel, colour = npar_phased>0)) + geom_point()
p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks()
p1

p1 <- ggplot(data = allsums, mapping = aes(x = tot_hetero, y = nparallel, colour = log2(nparallel/tot_hetero)>-8)) + geom_point()
p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks()
p1

p1 <- ggplot(data = allsums, mapping = aes(x = log2(nparallel/tot_hetero), fill = npar_phased>0)) + geom_histogram()
p1


write.table(x = allsums, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums[which(allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1), ], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(x = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", quote = F, col.names = T, row.names = F, sep = "\t")
