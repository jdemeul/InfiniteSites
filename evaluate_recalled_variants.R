### QC of PCAWG recalls compare to muts in consensus, report total cons, both (%), unique M2 (%)


library(VariantAnnotation)
# library(GenomicRanges)
library(parallel)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/code/check_germline_artefact.R") # for added QC
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R") # for added QC and release table

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
recalledsamples <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/samples_to_recall_success.txt", as.is = T)$tumor_wgs_aliquot_id
RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
BAMDIR <- "/srv/shared/vanloo/ICGC/"

releasetable <- read_pcawg_release_table(release_table_file = RELEASETABLEFILE)

# sampleid <- recalledsamples[which(recalledsamples == "0980e7fd-051d-45e9-9ca6-2baf073da4e8")]


get_caller_stats <- function(sampleid, releasetable) {
  
  print(paste0(sampleid))
  
  # read the new calls, subset to standard chroms, PASS variant and remove deletions
  vcffilepath <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/ICGC-infinite-sites-Mutect2Recallpanel/", pattern = paste0(sampleid, "_tumor_mutect2_snvs_indels_twicefiltered.vcf.gz$"), full.names = T, recursive = T)
  snvs <- readVcf(file = vcffilepath, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
  snvs <- snvs[seqnames(snvs) %in% c(1:22, "X", "Y") & rowRanges(snvs)$FILTER == "PASS" & nchar(ref(snvs)) == 1]
  
  # set aside biallelics
  snvs_biallelic <- snvs[which(lengths(alt(snvs)) == 2)]
  
  # remove insertions
  snvs <- snvs[lengths(alt(snvs)) == 1]
  snvs <- snvs[unlist(nchar(alt(snvs))) == 1]
  
  # read PCAWG calls
  snvs_pcawg <- readVcf(file = paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz"), genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))

  # only SNVs (i.e. clean up biallelics for single base subs, with no hint of germline)
  if (length(snvs_biallelic) > 0) {
    if (!is.null(info(snvs_biallelic)$P_GERMLINE)) {
      is_clean_biallelic <- which(sapply(X = alt(snvs_biallelic), FUN = function(x) all(x %in% c("A", "C", "G", "T"))) &
                                    sapply(X = info(snvs_biallelic)$P_GERMLINE, FUN = function(x) all(x < -1)) &
                                    sapply(X = info(snvs_biallelic)$NLOD, FUN = function(x) all(x > 2.2)))
    } else {
      is_clean_biallelic <- which(sapply(X = alt(snvs_biallelic), FUN = function(x) all(x %in% c("A", "C", "G", "T"))) &
                                    sapply(X = info(snvs_biallelic)$NLOD, FUN = function(x) all(x > 2.2)))
    }
    snvs_biallelic <- snvs_biallelic[is_clean_biallelic]
  }
  
  # turn into simple GRanges obj
  snvs_biallelic_temp <- rowRanges(snvs_biallelic)
  if (length(snvs_biallelic_temp) > 0) {
    tidx <- which(colnames(geno(snvs_biallelic)$F1R2) == sampleid)
    adddf <- data.frame(do.call(rbind, geno(snvs_biallelic)$F1R2[,-tidx]) + do.call(rbind, geno(snvs_biallelic)$F2R1[,-tidx]), 
                        do.call(rbind, geno(snvs_biallelic)$F1R2[,tidx]) + do.call(rbind, geno(snvs_biallelic)$F2R1[,tidx]))
  } else {
    adddf <- data.frame(cbind(integer(), integer(), integer(), integer(), integer(), integer()))
  }
  colnames(adddf) <- c(paste0("n_", c("ref", "alt1", "alt2"), "_reads"), paste0("t_", c("ref", "alt1", "alt2"), "_reads"))
  mcols(snvs_biallelic_temp) <- cbind(mcols(snvs_biallelic_temp), info(snvs_biallelic), adddf)
  snvs_biallelic <- snvs_biallelic_temp # overwrite
  rm(snvs_biallelic_temp)
  
  
  ### annotate variants further for QC with c("n_ref_count", "n_alt_count", "n_total_cov", "hg38clean", "Validation_status") # added 2021 
  if (length(snvs_biallelic) > 0) {
    snvs_biallelic <- add_snv_qc(sampleid = sampleid,
                                 sampledir = dirname(vcffilepath),
                                 releasetable = releasetable,
                                 snvs_vcf = snvs_biallelic,
                                 hg19tohg38 = hg19tohg38chain,
                                 ncores = 1,
                                 bamdir = "/srv/shared/vanloo/ICGC-infinite-sites-Mutect2Recallpanel/")
    snvs_biallelic <- snvs_biallelic[snvs_biallelic$n_total_cov >= 19 & snvs_biallelic$hg38clean]
  } else {
    mcols(snvs_biallelic)[, c("n_ref_count", "n_alt_count", "n_total_cov", "n_alt_count2", "Validation_status", "hg38clean")] <- 
      data.frame("n_alt_count" = integer(), "n_total_cov" = integer(), "n_alt_count2" = integer(), "Validation_status" = logical(), "hg38clean" = logical())
  }

  
  # annotate biallelics with allele-specific CN
  cn_pcawg <- GRanges(read.delim(file = paste0("/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/", sampleid, ".consensus.20170119.somatic.cna.annotated.txt"), as.is = T)[, 1:10], seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
  snvs_biallelic_df <- as.data.frame(snvs_biallelic)
  snvs_biallelic_df$alt1 <- as.character(DNAStringSet(sapply(X = snvs_biallelic_df$ALT, FUN = function(g) g[[1]])))
  snvs_biallelic_df$alt2 <- as.character(DNAStringSet(sapply(X = snvs_biallelic_df$ALT, FUN = function(g) g[[2]])))
  snvs_biallelic_df <- snvs_biallelic_df[, -8]
  snvs_biallelic_df[, c("total_cn", "major_cn", "minor_cn")] <- as.data.frame(mcols(cn_pcawg)[nearest(x = snvs_biallelic, subject = cn_pcawg), c("total_cn", "major_cn", "minor_cn")])
  
  # annotate which allele called in PCAWG
  snvs_biallelic_df$in_PCAWG <- rep("", nrow(snvs_biallelic_df))
  matchidxs <- match(x = paste0(snvs_biallelic_df$seqnames, "_", snvs_biallelic_df$start), table = paste0(seqnames(snvs_pcawg), "_", start(snvs_pcawg)))
  if (any(!is.na(matchidxs))) {
    snvs_biallelic_df[which(!is.na(matchidxs)), "in_PCAWG"] <- as.character(unlist(alt(snvs_pcawg))[matchidxs[!is.na(matchidxs)]])
  }

  # general stats
  # ninboth <- sum(overlapsAny(query = snvs, subject = snvs_pcawg, type = "equal")) + sum(snvs_biallelic_df$in_PCAWG != "")
  # outv <- c(onlym2 = length(snvs) - ninboth + 2*sum(snvs_biallelic_df$in_PCAWG == "") + sum(snvs_biallelic_df$in_PCAWG != ""),
  #           both = ninboth,
  #           onlycons = length(snvs_pcawg) - ninboth,
  #           nbiallelics = length(snvs_biallelic))

  ninboth <- sum(overlapsAny(query = snvs, subject = snvs_pcawg, type = "equal"))
  outv <- c(onlym2 = length(snvs) - ninboth + 2*sum(snvs_biallelic_df$in_PCAWG == "") + sum(snvs_biallelic_df$in_PCAWG != ""),
            both = ninboth + sum(snvs_biallelic_df$in_PCAWG != ""),
            onlycons = length(snvs_pcawg) - ninboth - sum(snvs_biallelic_df$in_PCAWG != ""),
            nbiallelics = nrow(snvs_biallelic_df))
  
  return(list(outv = outv, snvs_biallelic = snvs_biallelic_df))
}

# undebug(get_caller_stats)
# debug(get_caller_stats)
# debug(add_snv_qc)
# temp2 <- get_caller_stats(sampleid = "14c5b81d-da49-4db1-9834-77711c2b1d38", releasetable = releasetable)
# recalledsamples[sapply(callerstats_all, is.null)]
# sampleid <- "0dd0718d-5ddf-4c59-8c47-0f51303daeb5"

callerstats_all <- lapply(X = recalledsamples, FUN = get_caller_stats, releasetable = releasetable)
# callerstats_all <- mclapply(X = recalledsamples, FUN = get_caller_stats, mc.cores = 2, mc.preschedule = F, releasetable = releasetable)
saveRDS(object = callerstats_all, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_list.RDS")
# callerstats_all <- readRDS("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_list.RDS")
# callervariants <- do.call(rbind, lapply(callerstats_all, FUN = function(x) {
#   y <- as.data.frame(rowRanges(x$snvs_biallelic))
#   y$alt1 <- as.character(DNAStringSet(sapply(X = y$ALT, FUN = function(g) g[[1]])))
#   y$alt2 <- as.character(DNAStringSet(sapply(X = y$ALT, FUN = function(g) g[[2]])))
#   y <- y[, -8]
#   return(y)
#   }))

# FIX
callervariantsls <- lapply(callerstats_all, FUN = function(x) {
  x$snvs_biallelic
})

sapply(callervariantsls, FUN = function(x) length(colnames(x)))
jointcols <- intersect(x = colnames(callervariantsls[[1]]), y = colnames(callervariantsls[[2]]))
# sapply(X = callervariantsls, FUN = function(x, colsel) setdiff(y = colnames(x), x = colsel), colsel = jointcols)

callervariants <- do.call(rbind, lapply(callervariantsls, FUN = function(x, colsel) x[, colsel], colsel = jointcols))

callervariants$sampleid <- rep(recalledsamples, sapply(callervariantsls, FUN = nrow))

# some additonal checks 2021
sumtab <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)
callervariants$histology_abbreviation <- sumtab[match(x = callervariants$sampleid, table = sumtab$samplename), "histology_abbreviation"]
# callervariants$clean <- sapply(callervariants$NLOD, FUN = function(x) all(x > 2.2)) &
#   callervariants$n_total_cov >= 20 &
#   callervariants$hg38clean
# # hist(callervariants$n_total_cov)
p1 <- ggplot(data = callervariants) + geom_jitter(mapping = aes(x = histology_abbreviation, y = n_total_cov), alpha = .3)
p1 <- p1 + scale_y_log10() + annotation_logticks(sides = "l") + theme(axis.text.x = element_text(angle =90))
p1
# hist(rowSums(callervariants[callervariants$clean, c("n_alt_count", "n_alt_count2")]))
sum(callervariants$clean)


write.table(x = callervariants, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_20210217.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# rereadtab <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_20210217.txt", as.is = T)

callerstats <- data.frame(do.call(rbind, lapply(callerstats_all, FUN = function(x) x$outv)))
callerstats$sampleid <- recalledsamples
# callerstats_noartfilt <- callerstats

# go for twice filtered data (seq art filtered), actually affects only 20 samples, mostly small effects, few bigger
# allstats <- cbind(callerstats_noartfilt[,-4],callerstats)
# colnames(allstats) <- c("f1_onlym2", "f1_both", "f1_onlycons", "f2_onlym2", "f2_both", "f2_onlycons", "sampleid")
# 
# p1 <- ggplot(allstats, data = aes(x= 1:nrow(allstats), y = ))
# 
# which(allstats$f1_onlym2 != allstats$f2_onlym2)

callerstats$histology_abbreviation <- sumtab[match(x = callerstats$sampleid, table = sumtab$samplename), "histology_abbreviation"]

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(callerstats$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(callerstats$histology_abbreviation)
# cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Skin-Melanoma-Mucosal", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", "#000000", '#FF4500', '#FF4500')


callerstats[order(callerstats$both+callerstats$onlycons, decreasing = T), "rank"] <- 1:nrow(callerstats)

p1 <- ggplot(data = callerstats, mapping = aes(x = rank, colour = histology_abbreviation)) + geom_point(mapping = aes(y = both/(both+onlycons)), shape = 20, alpha = .8, show.legend = F) + 
  geom_point(mapping = aes(y = onlym2/(both+onlycons)), shape = 1, alpha = .8, show.legend = F)
p1 <- p1 + scale_color_manual(values = cvect) + ylim(c(0,1))
p1 <- p1 + theme_minimal() + theme(panel.grid.minor = element_blank()) + labs(y = "Fraction of consensus SNVs")
p1


p2 <- ggplot(data = callerstats, mapping = aes(x = rank, y = both+onlycons, colour = histology_abbreviation)) + geom_point(alpha = .8, show.legend = F, shape = 20) + scale_y_log10()
p2 <- p2 + scale_color_manual(values = cvect)
p2 <- p2 + theme_minimal() + theme(panel.grid.minor = element_blank()) + annotation_logticks(sides = "l", base = 10) + labs(y = "# consensus SNVs")
p2

ggsave(plot = p1, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Mutect2Valid_TP_20210217.pdf", width = 9, height = 2)
ggsave(plot = p2, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Mutect2Valid_ML_20210217.pdf", width = 9, height = 2)

### will want to exclude 142b6dbf-c943-4a7d-8ab6-13a975f48d7a due to high number of false positives (calls in clipped reads apparently)


p2 <- ggplot(data = callerstats, mapping = aes(x = rank, y = both+onlycons, colour = histology_abbreviation)) + geom_point(alpha = .8, show.legend = F, shape = 20) + scale_y_log10()
p2 <- p2 + geom_point(mapping = aes(y = nbiallelics), alpha = .8, show.legend = F, shape = 17) + scale_color_manual(values = cvect)
p2 <- p2 + theme_minimal() + theme(panel.grid.minor = element_blank()) + annotation_logticks(sides = "l", base = 10) + labs(y = "# consensus SNVs")
p2
ggsave(plot = p2, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Mutect2Valid_ML_plusThirds_20210217.pdf", width = 9, height = 2)


write.table(x = callerstats, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary_20210217.txt", quote = F, col.names = T, row.names = F, sep = "\t")


quantile(callerstats$both/(callerstats$both+callerstats$onlycons), probs = c(.025, .5, .975))
quantile(callerstats$onlym2/(callerstats$both+callerstats$onlycons), probs = c(.025, .5, .975))
quantile(2*callerstats$nbiallelics/(callerstats$onlym2+callerstats$both), probs = c(.025, .5, .975))
quantile(2*callerstats$nbiallelics/(callerstats$onlym2+callerstats$both), probs = c(.025, .5, .975))

nnewalleles <- c(by(data = callervariants$in_PCAWG, INDICES = callervariants$sampleid, FUN = function(x) sum(x != "")+ 2*sum(x == "")))
fracnewalleles <- nnewalleles[callerstats$sampleid]/callerstats$onlym2
fracnewalleles[is.na(fracnewalleles)] <- 0
quantile(fracnewalleles, probs = c(.025, .5, .975))*100

sum(callervariants$in_PCAWG != "")/nrow(callervariants)

sum(callervariants$in_PCAWG != "")+ 2*sum(callervariants$in_PCAWG == "")
