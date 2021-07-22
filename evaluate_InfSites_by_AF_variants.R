#### evaluate parallel calls

library(ggplot2)

source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
CLEANHISTOLOGYFILE <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

allsumfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210702_vafpipeline_out_alphapt1_hetonly/", pattern = "_InfSites_VAFpipeline_summarystats.txt", full.names = T, recursive = T)
allsums <- lapply(X = allsumfiles, FUN = read.table, row.names = 1, as.is = T)
allsums <- data.frame(t(do.call(cbind, allsums)))
allsums$sampleid <- gsub(pattern = "_InfSites_VAFpipeline_summarystats.txt", replacement = "", x = basename(allsumfiles))

allpcfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210702_vafpipeline_out_alphapt1_hetonly/", pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
allpc <- lapply(X = allpcfiles, FUN = read.table, col.names = c("pseudocount", "snv_slope"), as.is = T)
allpc <- data.frame(do.call(rbind, allpc))
allpc$sampleid <- gsub(pattern = "_pseudocount_calibrated.txt", replacement = "", x = basename(allpcfiles))

allsums <- merge(x = allsums, y = allpc, by = "sampleid")

sumtab <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)
allsums$histology_abbreviation <- sumtab[match(x = allsums$sampleid, table = sumtab$samplename), "histology_abbreviation"]
consrhopsi <- read.delim(file = "/camp/project/proj-emedlab-vanloo/ICGC_consensus_copynumber/consensus.20170217.purity.ploidy.txt.gz", as.is = T)
allsums <- cbind(allsums, consrhopsi[match(x = allsums$sampleid, table = consrhopsi$samplename), -1])
allsums$is_preferred <- allsums$sampleid %in% sumtab[sumtab$is_preferred, "samplename"]


allsimfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_all_writefrac/", pattern = "_infsites_permut_effgenfrac1.txt", full.names = T, recursive = T)
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




### note
RELEASETABLEFILE <- "/camp/project/proj-emedlab-vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
reltable <- read.delim(file = RELEASETABLEFILE, as.is = T)
sampleassignments <- strsplit(x = reltable$tumor_wgs_aliquot_id, split = ",")
sampleassignmentsv <- setNames(object = rep(x = reltable$wgs_exclusion_white_gray, lengths(sampleassignments)), nm = unlist(sampleassignments))
#note

resamplingmodel <- read.delim(file = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA_NODRIV.txt", as.is = T)
# resamplingmodel <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA.txt", as.is = T)
divergent_variants <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary_20210217.txt", as.is = T)


allsums$simnumber_new <- resamplingmodel[match(x = allsums$sampleid, table = resamplingmodel$sampleid), "total_new"]
allsums$ndivergent <- divergent_variants[match(x = allsums$sampleid, table = divergent_variants$sampleid), "nbiallelics"]
allsums$exclusionlist <- sampleassignmentsv[allsums$sampleid]



##### READING AND WRITING OUT THE HITS
write.table(x = allsums, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_20210701.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(x = allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope <= 1), ],
            file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210701.txt", quote = F, col.names = T, row.names = F, sep = "\t")

### writing out the specific variants for those samples which passed required filters
parhitsamples <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210701.txt", as.is = T)

# read ACTUAL hits
parhitsfiles <- paste0("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210702_vafpipeline_out_alphapt1_hetonly/", parhitsamples$sampleid, "/", parhitsamples$sampleid, "_snv_mnv_infSites_finalhits.txt")
parhits <- lapply(X = parhitsfiles, FUN = read.delim, as.is = T)
parhitsdf <- do.call(rbind, parhits)
parhitsdf$sampleid <- rep(parhitsamples$sampleid, sapply(parhits, nrow))


parhitsdf$ref <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$ref))
parhitsdf$alt <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$alt))

write.table(x = parhitsdf, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210701_finalset.txt", quote = F, row.names = F, col.names = T, sep = "\t")



allsums <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_20210701.txt", as.is = T)

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(allsums$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(allsums$histology_abbreviation)


# clean set of samples for PLOTTING too
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
# ggsave(plot = p1, filename = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedPermutation.pdf", width = 15, height = 12)
ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedPermutation.pdf", width = 15, height = 12)




p1 <- ggplot(data = df1, mapping = aes(x = simnumber_new/1000, y = nparallel+ndivergent))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_smooth() + geom_point(shape = 21, mapping = aes(fill = histology_abbreviation), color = "black", alpha = .6, size = 5) + 
  annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_newsim$estimate[[1]], digits = 3), collapse = "")) +
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3)) + annotation_logticks() + scale_fill_manual(values = cvect, name = "Cancer type") + 
  labs(x = "Average biallelic violations per 1000 simulations (Resampling model)", y = "Number of biallelic variants detected (divergent + parallel)") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4), ylim = c(1, 1e4))
p1
ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedResampling_NODRIVERS.pdf", width = 15, height = 12)
write.table(x = df1, file = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulated.txt_NODRIVERS", quote = F, row.names = F, sep = "\t")


### added one to use the estimate based on phasing to the plot
cv_newsim_mod <- cor.test(df1$simnumber_new, ifelse(df1$tot_phaseable>1e4, df1$med, df1$nparallel) + df1$ndivergent)

p1 <- ggplot(data = df1, mapping = aes(x = simnumber_new/1000, y = ifelse(tot_phaseable>1e4, med, nparallel) +ndivergent))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_smooth() + geom_point(shape = 21, mapping = aes(fill = histology_abbreviation, color = tot_phaseable>1e4, size = tot_hetero), alpha = .6, size = 5) + 
  scale_color_manual(values = c('TRUE' = "red", 'FALSE' = "black"), name = "# phaseable SNVs > 10,000") + scale_size()
  annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_newsim_mod$estimate[[1]], digits = 3), collapse = "")) +
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3)) + annotation_logticks() + scale_fill_manual(values = cvect, name = "Cancer type") + 
  labs(x = "Average biallelic violations per 1000 simulations (Resampling model)", y = "Number of biallelic variants detected/estimated (divergent + parallel)") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4), ylim = c(1, 1e4))
# p1 <- p1 + geom_point(data = df1[df1$tot_phaseable>1e4,], mapping = aes(x = simnumber_new/1000, y = med+ndivergent, fill = histology_abbreviation), alpha = .1)
p1
# ggsave(plot = p1, filename = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/ObservedOrEstimated_biallelic_variants_vs_simulatedResampling.pdf", width = 15, height = 12)
ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/ObservedOrEstimated_biallelic_variants_vs_simulatedResampling_NODRIVER.pdf", width = 15, height = 12)


