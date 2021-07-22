#check parsite muatbility #3


# load all variants
library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
# library(ggplot2)
library(VariantAnnotation)
library(parallel)
# library(rslurm)
# library(VGAM)
library(ggplot2)
library(ggrepel)

# source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")
# source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")

SNVMNVINDELDIR <- "/camp/project/proj-vanloo-secure/PCAWG/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/"
# RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
# sampleid <- "fc876f5c-8339-bc9c-e040-11ac0d485160"
# sampleid <- "7cae6c0b-36fe-411b-bbba-093a4c846d84"
# sampleid <- "0332b017-17d5-4083-8fc4-9d6f8fdbbbde"
# sampleid <- "804ffa2e-158b-447d-945c-707684134c87"

sumtab <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
sampleids <- sumtab[sumtab$is_preferred, "samplename"]

parhits <- GRanges(read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210701_finalset.txt", as.is = T))
parhits <- parhits[!seqnames(parhits) == "Y"]
parhits <- parhits[parhits$sampleid %in% sampleids]
thirdhits <- GRanges(read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T))
thirdhits <- thirdhits[!seqnames(thirdhits) == "Y"]
thirdhits <- thirdhits[thirdhits$sampleid %in% sampleids]

# only add those loci which are NOT already in there (with same sample)
thirdhitsforaddition <- granges(thirdhits, use.names = F, use.mcols = F)[thirdhits$in_PCAWG == ""]
mcols(thirdhitsforaddition) <- DataFrame(trinuc = getSeq(resize(thirdhitsforaddition, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5))
# thirdhitsforaddition <- c(granges(thirdhits),granges(thirdhits))
# mcols(thirdhitsforaddition) <- DataFrame(alt = c(thirdhits$alt1,thirdhits$alt2),
#                                          trinuc = rep(getSeq(resize(thirdhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5), 2),
#                                          ref = rep(thirdhits$REF, 2))
flipthirdidx <- which(subseq(thirdhitsforaddition$trinuc, 2,2) %in% c("A", "G"))
thirdhitsforaddition$trinuc[flipthirdidx] <- reverseComplement(thirdhitsforaddition$trinuc[flipthirdidx])
# thirdhitsforaddition$alt <- DNAStringSet(thirdhitsforaddition$alt)
# thirdhitsforaddition$alt[flipthirdidx] <- reverseComplement(thirdhitsforaddition$alt[flipthirdidx])


driverhits <- read.delim(file = "/camp/project/proj-vanloo-secure/PCAWG/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
driverhits <- driverhits[driverhits$ref %in% c("A", "C", "G", "T") & driverhits$alt %in% c("A", "C", "G", "T"), c("chr", "pos")]
driverhits <- sort(unique(GRanges(seqnames = driverhits$chr, ranges = IRanges(start = driverhits$pos, width = 1))))

# sum(sapply(sampleids, FUN = function(x, snvsmnvdir) file.exists(file.path(snvsmnvdir, pattern = paste0(x, ".consensus.20160830.somatic.snv_mnv.vcf.gz"))), snvsmnvdir = SNVMNVINDELDIR))
# 2583/2658 samples

get_snvs <- function(sampleid, snvsmnvdir) {
  print(paste0("Loading snvs for sample ", sampleid))
  snv_mnvfile <- file.path(snvsmnvdir, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz"))
  
  if (!file.exists(snv_mnvfile)) {
    snvs <- GRanges()
    snvs$alt <- character()
  } else {
    snvs <- readVcf(file = snv_mnvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
    snvs <- granges(rowRanges(snvs[lengths(alt(snvs)) == 1 & lengths(ref(snvs)) == 1 & seqnames(snvs) != "Y"]), use.names = F, use.mcols = F)
    # mcols(snvs) <- DataFrame(alt = as.character(unlist(snvs$ALT)))
  }
  
  snvs$trinuc <- getSeq(resize(snvs, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
  # snvs$ref <- subseq(snvs$trinuc, 2,2)
  flipidx <- which(subseq(snvs$trinuc, 2,2) %in% c("A", "G"))
  snvs$trinuc[flipidx] <- reverseComplement(snvs$trinuc[flipidx])
  # snvs$alt <- DNAStringSet(snvs$alt)
  # snvs$alt[flipidx] <- reverseComplement(snvs$alt[flipidx])
  
  # snvs <- snvs[which(snvs$trinuc == "TCT")]
  
  return(snvs)
}

### remove all driver muts to reduce influence of selection
# debug(get_snvs)
allsnvs <- unlist(GRangesList(mclapply(X = sampleids, FUN = get_snvs, snvsmnvdir = SNVMNVINDELDIR, mc.preschedule = T, mc.cores = 6)))
allsnvs <- c(allsnvs, thirdhitsforaddition)#[which(thirdhitsforaddition$trinuc == "TCT")])
# allsnvs$trinucmut <- factor(paste0(allsnvs$trinuc, ">", allsnvs$alt))
allsnvs <- subsetByOverlaps(x = allsnvs, ranges = driverhits, invert = T)

saveRDS(object = allsnvs, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/mutability_check_allsnvs+bial_20210701.rds")
# allsnvs <- readRDS(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/mutability_check_allsnvs+bial.rds")

# allsnvs$trinuc <- getSeq(resize(allsnvs, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
# allsnvs$ref <- subseq(allsnvs$trinuc, 2,2)
# flipidx <- which(allsnvs$ref %in% c("A", "G"))
# allsnvs$trinuc[flipidx] <- reverseComplement(allsnvs$trinuc[flipidx])
# allsnvs$alt <- DNAStringSet(allsnvs$alt)
# allsnvs$alt[flipidx] <- reverseComplement(allsnvs$alt[flipidx])
# allsnvs$trinucmut <- factor(paste0(allsnvs$trinuc, ">", allsnvs$alt))

# do per trinuc type
allsnvs_dedup <- granges(x = sort(unique(allsnvs)), use.mcols = T, use.names = F)
# mcols(allsnvs_dedup) <- DataFrame(trinuc = allsnvs_dedup$trinuc)
# colnames(mcols(allsnvs_dedup)) <- "trinuc"

# allsnvs <- granges(allsnvs, use.names=FALSE, use.mcols=FALSE)

allsnvs_dedup$is_parhit <- as.factor(overlapsAny(query = allsnvs_dedup, subject = c(granges(parhits), granges(thirdhits)), type = "equal"))
allsnvs_dedup$hitcount <- countOverlaps(query = allsnvs_dedup, subject = allsnvs, type = "equal")
allsnvs_dedup$trinuc <- as.factor(allsnvs_dedup$trinuc)
allsnvs_dedup_dfonly <- mcols(allsnvs_dedup, use.names = F)[, c("trinuc", "is_parhit", "hitcount")]

rm(allsnvs, allsnvs_dedup)
gc()

get_bial_freq_strat <- function(alloci) {
  out <- aggregate(x = alloci$is_parhit,
                   by = list(trinuc = alloci$trinuc, hitcount = alloci$hitcount),
                   FUN = table, simplify = F, drop = T)
  outdf <- data.frame(out[, 1:2], do.call(rbind, out[,3]))
  return(outdf)
}

# debug(get_bial_freq_strat)
outdf <- get_bial_freq_strat(alloci = allsnvs_dedup_dfonly)
matchidx <- paste0(outdf$trinuc, "_", outdf$hitcount)

# outdf[outdf$trinuc == "TCT", ]

# partrinucfreqs <- sort(table(allsnvs_dedup$trinuc[which(allsnvs_dedup$is_parhit)]))
# selecttrinucs <- names(which(partrinucfreqs >= 100))

# RNGkind("L'Ecuyer-CMRG")
# set.seed(953755)
# mc.reset.stream()

replicdfs <- mclapply(X = 1:100, FUN = function(x) get_bial_freq_strat(alloci = allsnvs_dedup_dfonly[sample(1:nrow(allsnvs_dedup_dfonly), replace = T), ]), mc.preschedule = T, mc.cores = 4)

outdf <- cbind(outdf, t(apply(X = do.call(cbind, lapply(X = replicdfs, FUN = function(x, idxv) {
  y <- x[match(x = idxv, table = paste0(x[,1], "_", x[,2])), 3:4]
  z <- y[, 2]/rowSums(y)
  z[is.na(z)] <- 0
  return(z)
  }, idxv = matchidx)), MARGIN = 1, FUN = quantile, probs = c(0.025, 0.5, 0.975))))

colnames(outdf) <- c("trinuc", "hitcount", "nr_mono", "nr_bi", "ci_lo", "ci_med", "ci_hi")

saveRDS(object = replicdfs, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/mutability_check_bootstraps_20210701.rds")
# replicdfs <- readRDS(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/mutability_check_bootstraps.rds")
# mclapply(X = 1:16, FUN = function(x) sample(1:10, replace = T), mc.preschedule = T, mc.cores = 16)

# singleloci <- as.data.frame(t(apply(X = do.call(cbind, lapply(X = replicdfs, FUN = function(x) x[, 4]/rowSums(x[, 3:4]))), MARGIN = 1, FUN = quantile, probs = c(0.025, 0.5, 0.975), na.rm = T)))
# replicdf <- cbind(replicdfs[[1]][, 1:2], singleloci)

# replicdf2 <- merge(x = replicdf, y = outdf, all = T)

write.table(x = outdf, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/mutation_rate_vs_mutability_bootstrapres_20210701.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# outdf <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/mutation_rate_vs_mutability_bootstrapres.txt", as.is = T)

plotdf <- outdf[outdf$nr_bi >= 5 & outdf$trinuc %in% outdf$trinuc[outdf$hitcount == 2 & outdf$nr_bi >= 5], ]
textdf <- plotdf[!duplicated(plotdf$trinuc, fromLast = T), ]

p1 <- ggplot(data = plotdf, mapping = aes(x = hitcount, colour = trinuc, fill = trinuc, y = nr_bi/(nr_mono+nr_bi))) + geom_path(alpha = .4) + geom_point()# + geom_pointrange(mapping = aes(ymin = ci_lo, ymax = ci_hi))
p1 <- p1 + geom_ribbon(mapping = aes(ymin = ci_lo, ymax = ci_hi), colour = NA, alpha = .1) + scale_y_log10() + annotation_logticks(sides = "l") + coord_cartesian(xlim = c(0.75,7))
p1 <- p1 + geom_text_repel(data = textdf, mapping = aes(label = trinuc), alpha = .75, nudge_x = .25, box.padding = 0, point.padding = 0, segment.colour = NA)
p1 <- p1 + theme_minimal() + theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "none") + scale_x_continuous(breaks = 1:7)
p1 <- p1 + labs(x = "# SNVs at loci", y = "Frequency of biallelic hits at loci")
p1

ggsave(filename = paste0("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/mutation_rate_vs_biallelic_rate_20210701.pdf"), plot = p1, width = 10, height = 5, useDingbats=FALSE)






