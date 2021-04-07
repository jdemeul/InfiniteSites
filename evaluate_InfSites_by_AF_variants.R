#### evaluate parallel calls

library(ggplot2)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

# allsumfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", pattern = "_InfSites_VAFpipeline_summarystats.txt", full.names = T, recursive = T)
allsumfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", pattern = "_InfSites_VAFpipeline_summarystats.txt", full.names = T, recursive = T)
allsums <- lapply(X = allsumfiles, FUN = read.table, row.names = 1, as.is = T)
allsums <- data.frame(t(do.call(cbind, allsums)))
allsums$sampleid <- gsub(pattern = "_InfSites_VAFpipeline_summarystats.txt", replacement = "", x = basename(allsumfiles))

# allpcfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
allpcfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly//", pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
allpc <- lapply(X = allpcfiles, FUN = read.table, col.names = c("pseudocount", "snv_slope"), as.is = T)
allpc <- data.frame(do.call(rbind, allpc))
allpc$sampleid <- gsub(pattern = "_pseudocount_calibrated.txt", replacement = "", x = basename(allpcfiles))

allsums <- merge(x = allsums, y = allpc, by = "sampleid")

sumtab <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)
allsums$histology_abbreviation <- sumtab[match(x = allsums$sampleid, table = sumtab$samplename), "histology_abbreviation"]
consrhopsi <- read.delim(file = "/srv/shared/vanloo/ICGC_consensus_copynumber/consensus.20170217.purity.ploidy.txt.gz", as.is = T)
allsums <- cbind(allsums, consrhopsi[match(x = allsums$sampleid, table = consrhopsi$samplename), -1])
allsums$is_preferred <- allsums$sampleid %in% sumtab[sumtab$is_preferred, "samplename"]


allsimfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_all_writefrac/", pattern = "_infsites_permut_effgenfrac1.txt", full.names = T, recursive = T)
allsims <- lapply(X = allsimfiles, FUN = function(x) {
  y <- read.table(x, as.is = T, header = T)
  nparsim <- sum(y[y$type == "parallel", "counts_tot"])/1000
  ndivsim <- sum(y[y$type == "third_allele", "counts_tot"])/1000
  return(c(nparsim = nparsim, ndivsim = ndivsim))
  })
allsims <- data.frame(do.call(rbind, allsims))
allsims$sampleid <- gsub(pattern = "_infsites_permut_effgenfrac1.txt", replacement = "", x = basename(allsimfiles))

allsums$simnumber <- allsims[match(x = allsums$sampleid, table = allsims$sampleid), "nparsim"]
allsums$simnumber_divergent <- allsims[match(x = allsums$sampleid, table = allsims$sampleid), "ndivsim"]
allsums$snvs_total <- rowSums(sumtab[match(x = allsums$sampleid, table = sumtab$samplename), c("num_clonal", "num_subclonal")])

# allsums$hitratio <- log10(allsums$nparallel/allsums$tot_testable)
# allsums$hitratio2 <- log10(allsums$med/allsums$tot_testable)
# allsums$hitratio3 <- log10(allsums$tot_testable^2/allsums$nparallel)
# allsums$hitratio4 <- log10(allsums$tot_testable^2/allsums$med)
# allsums$hitratio5 <- log10(allsums$tot_diploid^2/allsums$med_diploid)
# allsums$pratio <- pnorm(q = allsums$hitratio3, mean = mean(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), sd = sd(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), lower.tail = TRUE, log.p = FALSE)
# allsums$simratio <- log2((allsums$nparallel+0.001)/(allsums$simnumber+0.001))
# allsums$simratio2 <- log2((allsums$nparallel+0.001)/(allsums$simnumber*allsums$tot_testable/allsums$snvs_total+0.001))

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(allsums$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(allsums$histology_abbreviation)


# temp <- read.delim("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_het_only_writefrac/0009b464-b376-4fbc-8a56-da538269a02f/0009b464-b376-4fbc-8a56-da538269a02f_infsites_permut_effgenfrac0.1.txt", as.is = T)
# temp <- sum(temp[temp$type == "parallel", "counts_tot"])/1000
# p1 <- ggplot(data = allsums[which(allsums$prec > 0 & allsums$rec > 0), ]) + geom_point(mapping = aes(x = hitratio2, y = purity_conf_mad))
# p1

# p1 <- ggplot(data = allsums) + geom_point(mapping = aes(x = hitratio, y = purity_conf_mad, colour = hitratio > -2.352241))
# p1

# ratecutoff <- qnorm(p = .01, mean = mean(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), sd = sd(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), lower.tail = TRUE, log.p = FALSE)
# p1 <- ggplot(data = allsums) + geom_point(mapping = aes(y = log10(nparallel), x = log10(snvs_total), color = histology_abbreviation, shape = log2((nparallel+0.001)/(simnumber+0.001)) < 8), size = 5) +
#   annotation_logticks() + scale_color_manual(values = cvect)
# p1


### note
RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
reltable <- read.delim(file = RELEASETABLEFILE, as.is = T)
sampleassignments <- strsplit(x = reltable$tumor_wgs_aliquot_id, split = ",")
sampleassignmentsv <- setNames(object = rep(x = reltable$wgs_exclusion_white_gray, lengths(sampleassignments)), nm = unlist(sampleassignments))
#note

resamplingmodel <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA.txt", as.is = T)
divergent_variants <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary_20210217.txt", as.is = T)


allsums$simnumber_new <- resamplingmodel[match(x = allsums$sampleid, table = resamplingmodel$sampleid), "total_new"]
allsums$ndivergent <- divergent_variants[match(x = allsums$sampleid, table = divergent_variants$sampleid), "nbiallelics"]
allsums$exclusionlist <- sampleassignmentsv[allsums$sampleid]

df1 <- allsums[which(allsums$snv_slope <= 1 & allsums$is_preferred), ]
# View(df1[is.na(df1$simnumber_new),])
# df1[is.na(df1$simnumber_new), "simnumber_new"] <- 0
df1[is.na(df1$ndivergent), "ndivergent"] <- 0

df1 <- df1[which(df1$simnumber+df1$simnumber_divergent >= 0.001 & (df1$simnumber_new >= 1 | is.na(df1$simnumber_new))), ]
cv_newsim <- cor.test(df1$simnumber_new, df1$nparallel+df1$ndivergent)
cv_oldsim <- cor.test(df1$simnumber+df1$simnumber_divergent, df1$nparallel+df1$ndivergent)

p1 <- ggplot(data = df1, mapping = aes(x = simnumber+simnumber_divergent, y = nparallel+ndivergent))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_smooth() + 
  geom_point(shape = 21, mapping = aes(fill = histology_abbreviation), color = "black", alpha = .6, size = 5) +
  # geom_point(shape = 21, mapping = aes(fill = exclusionlist), alpha = .6, size = 5) + 
  annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_oldsim$estimate[[1]], digits = 3), collapse = "")) +
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3)) + annotation_logticks() +
  scale_fill_manual(values = cvect, name = "Cancer type") +
  # geom_errorbar(mapping = aes(ymin = lower, ymax = upper)) +
  labs(x = "Average biallelic violations per 1000 simulations (Permutation model)", y = "Number of biallelic variants detected (divergent + parallel)") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4), ylim = c(1, 1e4))
p1
ggsave(plot = p1, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedPermutation.pdf", width = 15, height = 12)

p1 <- ggplot(data = df1, mapping = aes(x = simnumber_new/1000, y = nparallel+ndivergent))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_smooth() + geom_point(shape = 21, mapping = aes(fill = histology_abbreviation), color = "black", alpha = .6, size = 5) + 
  annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_newsim$estimate[[1]], digits = 3), collapse = "")) +
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3)) + annotation_logticks() + scale_fill_manual(values = cvect, name = "Cancer type") + 
  labs(x = "Average biallelic violations per 1000 simulations (Resampling model)", y = "Number of biallelic variants detected (divergent + parallel)") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4), ylim = c(1, 1e4))
p1
ggsave(plot = p1, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedResampling.pdf", width = 15, height = 12)
write.table(x = df1, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulated.txt", quote = F, row.names = F, sep = "\t")




# 

# # p1 <- ggplot(data = allsums) + geom_point(mapping = aes(y = simratio, x = simratio2))
# # p1
# 
# sum(log2((allsums$nparallel+0.001)/(allsums$simnumber+0.001)) > 8)
# sum(allsums$pratio < 0.01)
# 
# p1 <- ggplot(data = allsums) + geom_point(mapping = aes(y = log10(nparallel), x = log10(snvs_total), colour = pratio < 0.01))
# p1
# 
# p1 <- ggplot(data = allsums) + geom_point(mapping = aes(y = log10(nparallel), x = log10(simnumber)))
# p1
# 
# p1 <- ggplot(data = allsums, mapping = aes(y = log10(snvs_tot), x = log10(tot_testable))) + geom_point() + geom_smooth()
# p1
# 
# p1 <- ggplot(data = allsums) + geom_histogram(mapping = aes(x = hitratio3, fill = !is.na(allsums$prec) & allsums$prec > 0.5 & allsums$rec > 0))
# p1
# 
# p1 <- ggplot() + geom_histogram(data = allsums[allsums$tot_phaseable >= 1000, ], mapping = aes(x = hitratio4), fill = "blue", alpha = .5) +
#   geom_histogram(data = allsums, mapping = aes(x = hitratio3), fill = "red", alpha = .5)
# p1
# 
# p1 <- ggplot() + geom_histogram(data = allsums[allsums$tot_phaseable >= 1000, ], mapping = aes(x = hitratio4), fill = "blue", alpha = .5) +
#   geom_histogram(data = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], mapping = aes(x = hitratio3), fill = "red", alpha = .5)
# p1


# p1 <- ggplot() + geom_point(data = allsums[allsums$tot_phaseable > 1000, ], mapping = aes(y = log10(med), x = log10(tot_testable^2)), colour = "black") +
#   geom_point(data = allsums, mapping = aes(y = log10(nparallel), x = log10(tot_testable^2), colour = hitratio > -2.3))
# p1



#this is THE set to work with:
# View(allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ]) #ALT
# View(allsums[which(allsums$simratio < 8 & allsums$nparallel > 0 & allsums$snv_slope < 1), ]) #ALT
# View(allsums[which(allsums$nparallel > 0 & allsums$snv_slope < 1 & allsums$is_preferred), ]) #ALT
View(allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope < 1), ])
nrow(allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope < 1), ])
# 
# 
# p1 <- ggplot(data = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], mapping = aes(x = med, y = nparallel, colour = histology_abbreviation)) + geom_point()
# p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks() + geom_abline()
# p1
# 
# p1 <- ggplot(data = allsums[which(allsums$simratio < 8 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], mapping = aes(x = med, y = nparallel, colour = histology_abbreviation)) + geom_point()
# p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks() + geom_abline()
# p1
# 
# p1 <- ggplot(data = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], mapping = aes(x = tot_hetero, y = nparallel, colour = npar_phased>0)) + geom_point()
# p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks()
# p1
# 
# p1 <- ggplot(data = allsums, mapping = aes(x = tot_hetero, y = nparallel, colour = log2(nparallel/tot_hetero)>-8)) + geom_point()
# p1 <- p1 + scale_x_log10() + scale_y_log10() + annotation_logticks()
# p1
# 
# p1 <- ggplot(data = allsums, mapping = aes(x = log2(nparallel/tot_hetero), fill = npar_phased>0)) + geom_histogram()
# p1


write.table(x = allsums, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_2021.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums[which(allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1), ], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(x = allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope < 1), ],
            file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_2021.txt", quote = F, col.names = T, row.names = F, sep = "\t")


### writing out the specific variants
parhitsamples <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_2021.txt", as.is = T)
# parhitsamples <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)


# read actual hits
# parhitsfiles <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", parhitsamples$sampleid, "/", parhitsamples$sampleid, "_snv_mnv_infSites_finalhits.txt")
parhitsfiles <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", parhitsamples$sampleid, "/", parhitsamples$sampleid, "_snv_mnv_infSites_finalhits.txt")
parhits <- lapply(X = parhitsfiles, FUN = read.delim, as.is = T)
parhitsdf <- do.call(rbind, parhits)
parhitsdf$sampleid <- rep(parhitsamples$sampleid, sapply(parhits, nrow))


parhitsdf$ref <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$ref))
parhitsdf$alt <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$alt))

write.table(x = parhitsdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_2021.txt", quote = F, row.names = F, col.names = T, sep = "\t")




### checks comparing 2021 vs 2020 calls
head(parhitsamples)
parhits_old <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)
inold <- setdiff(x = parhits_old$sampleid, y = parhitsamples$sampleid)
innew <- setdiff(y = parhits_old$sampleid, x = parhitsamples$sampleid)

parhits_old[parhits_old$sampleid %in% inold, "nparallel"]
parhitsamples[parhitsamples$sampleid %in% innew, "nparallel"]

mergedf <- merge(x = parhits_old, y = parhitsamples, by = "sampleid", all = T, suffixes = c("_old", "_new"))

p1 <- ggplot(data = mergedf, mapping = aes(x = log10(nparallel_old), y = log10(nparallel_new))) + geom_point() + geom_abline() + annotation_logticks()
p1

p1 <- ggplot(data = mergedf, mapping = aes(x = prec_new-prec_old)) + geom_histogram()
p1
p1 <- ggplot(data = mergedf, mapping = aes(x = rec_new-rec_old)) + geom_histogram()
p1

sum(parhits_old$nparallel)
sum(parhitsamples$nparallel)

missingnow <- sumtab[!sumtab$samplename %in% allsums$sampleid, ]
sort(rowSums(sumtab[!sumtab$samplename %in% allsums$sampleid, c("num_clonal", "num_subclonal")]))

mergedf[which(mergedf$nparallel_dipl_old - mergedf$nparallel_dipl_new > 100), c("sampleid", "snvs_total", "nparallel_old", "nparallel_new", "prec_old", "prec_new", "rec_old", "rec_new")]

