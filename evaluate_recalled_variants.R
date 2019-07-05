### QC of PCAWG recalls compare to muts in consensus, report total cons, both (%), unique M2 (%)


library(VariantAnnotation)
# library(GenomicRanges)
library(parallel)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

recalledsamples <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/samples_to_recall_success.txt", as.is = T)$tumor_wgs_aliquot_id

# sampleid <- recalledsamples[which(recalledsamples == "0980e7fd-051d-45e9-9ca6-2baf073da4e8")]


get_caller_stats <- function(sampleid) {
  
  print(paste0(sampleid))
  
  vcffilepath <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/", pattern = paste0(sampleid, "_tumor_mutect2_snvs_indels_twicefiltered.vcf.gz$"), full.names = T, recursive = T)
  snvs <- readVcf(file = vcffilepath, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
  snvs <- snvs[seqnames(snvs) %in% c(1:22, "X", "Y") & rowRanges(snvs)$FILTER == "PASS" & nchar(ref(snvs)) == 1]
  
  # set aside biallelics
  snvs_biallelic <- snvs[which(lengths(alt(snvs)) == 2)]
  
  snvs <- snvs[lengths(alt(snvs)) == 1]
  snvs <- snvs[unlist(nchar(alt(snvs))) == 1]
  
  snvs_pcawg <- readVcf(file = paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz"), genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))

  # only SNVs
  if (length(snvs_biallelic) > 0) {
    if (!is.null(info(snvs_biallelic)$P_GERMLINE)) {
      is_clean_biallelic <- which(sapply(X = alt(snvs_biallelic), FUN = function(x) all(x %in% c("A", "C", "G", "T"))) &
                                    sapply(X = info(snvs_biallelic)$P_GERMLINE, FUN = function(x) all(x < -1)))
    } else {
      is_clean_biallelic <- which(sapply(X = alt(snvs_biallelic), FUN = function(x) all(x %in% c("A", "C", "G", "T"))))
    }
    snvs_biallelic <- snvs_biallelic[is_clean_biallelic]
  }
  
  # annotate biallelics with allele-specific CN
  cn_pcawg <- GRanges(read.delim(file = paste0("/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/", sampleid, ".consensus.20170119.somatic.cna.annotated.txt"), as.is = T)[, 1:10], seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
  snvs_biallelic_df <- as.data.frame(rowRanges(snvs_biallelic))
  snvs_biallelic_df$alt1 <- as.character(DNAStringSet(sapply(X = snvs_biallelic_df$ALT, FUN = function(g) g[[1]])))
  snvs_biallelic_df$alt2 <- as.character(DNAStringSet(sapply(X = snvs_biallelic_df$ALT, FUN = function(g) g[[2]])))
  snvs_biallelic_df <- snvs_biallelic_df[, -8]
  snvs_biallelic_df[, c("total_cn", "major_cn", "minor_cn")] <- as.data.frame(mcols(cn_pcawg)[nearest(x = snvs_biallelic, subject = cn_pcawg), c("total_cn", "major_cn", "minor_cn")])
  
  # general stats
  ninboth <- sum(overlapsAny(query = snvs, subject = snvs_pcawg, type = "equal"))
  outv <- c(onlym2 = length(snvs) - ninboth,
            both = ninboth,
            onlycons = length(snvs_pcawg) - ninboth,
            nbiallelics = length(snvs_biallelic))

  return(list(outv = outv, snvs_biallelic = snvs_biallelic_df))
}

# debug(get_caller_stats)
# callerstats_all <- lapply(X = recalledsamples[1:4], FUN = get_caller_stats)
callerstats_all <- mclapply(X = recalledsamples, FUN = get_caller_stats, mc.cores = 20, mc.preschedule = F)

# callervariants <- do.call(rbind, lapply(callerstats_all, FUN = function(x) {
#   y <- as.data.frame(rowRanges(x$snvs_biallelic))
#   y$alt1 <- as.character(DNAStringSet(sapply(X = y$ALT, FUN = function(g) g[[1]])))
#   y$alt2 <- as.character(DNAStringSet(sapply(X = y$ALT, FUN = function(g) g[[2]])))
#   y <- y[, -8]
#   return(y)
#   }))

callervariants <- do.call(rbind, lapply(callerstats_all, FUN = function(x) {
  x$snvs_biallelic
}))

callervariants$sampleid <- rep(recalledsamples, sapply(callerstats_all, FUN = function(x) nrow(x$snvs_biallelic)))
write.table(x = callervariants, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", sep = "\t", quote = F, col.names = T, row.names = F)


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

sumtab <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)
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

ggsave(plot = p1, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Mutect2Valid_TP.pdf", width = 9, height = 2)
ggsave(plot = p2, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Mutect2Valid_ML.pdf", width = 9, height = 2)

### will want to exclude 142b6dbf-c943-4a7d-8ab6-13a975f48d7a due to high number of false positives (calls in clipped reads apparently)


p2 <- ggplot(data = callerstats, mapping = aes(x = rank, y = both+onlycons, colour = histology_abbreviation)) + geom_point(alpha = .8, show.legend = F, shape = 20) + scale_y_log10()
p2 <- p2 + geom_point(mapping = aes(y = nbiallelics), alpha = .8, show.legend = F, shape = 17) + scale_color_manual(values = cvect)
p2 <- p2 + theme_minimal() + theme(panel.grid.minor = element_blank()) + annotation_logticks(sides = "l", base = 10) + labs(y = "# consensus SNVs")
p2
ggsave(plot = p2, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Mutect2Valid_ML_plusThirds.pdf", width = 9, height = 2)


write.table(x = callerstats, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary.txt", quote = F, col.names = T, row.names = F)

