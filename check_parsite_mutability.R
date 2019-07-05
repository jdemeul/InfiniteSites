### check mutation rates at sites with violations

# load all variants
library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
# library(ggplot2)
library(VariantAnnotation)
library(parallel)
# library(rslurm)
# library(VGAM)
# library(ggplot2)

# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")
# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")

SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/"
# RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
# sampleid <- "fc876f5c-8339-bc9c-e040-11ac0d485160"
# sampleid <- "7cae6c0b-36fe-411b-bbba-093a4c846d84"
# sampleid <- "0332b017-17d5-4083-8fc4-9d6f8fdbbbde"
# sampleid <- "804ffa2e-158b-447d-945c-707684134c87"

sumtab <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
sampleids <- sumtab[sumtab$is_preferred, "samplename"]

parhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt", as.is = T))
parhits <- parhits[!seqnames(parhits) == "Y"]
thirdhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T))
thirdhits <- thirdhits[!seqnames(thirdhits) == "Y"]

thirdhitsforaddition <- granges(thirdhits)
mcols(thirdhitsforaddition) <- DataFrame(alt = thirdhits$alt1,
                                         trinuc = getSeq(resize(thirdhitsforaddition, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5),
                                         ref = thirdhits$REF)
flipthirdidx <- which(thirdhitsforaddition$ref %in% c("A", "G"))
thirdhitsforaddition$trinuc[flipthirdidx] <- reverseComplement(thirdhitsforaddition$trinuc[flipthirdidx])
thirdhitsforaddition$alt <- DNAStringSet(thirdhitsforaddition$alt)
thirdhitsforaddition$alt[flipthirdidx] <- reverseComplement(thirdhitsforaddition$alt[flipthirdidx])


driverhits <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
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
    snvs <- rowRanges(readVcf(file = snv_mnvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
    snvs <- snvs[lengths(mcols(snvs)$ALT) == 1 & lengths(mcols(snvs)$REF) == 1 & seqnames(snvs) != "Y"]
    mcols(snvs) <- DataFrame(alt = as.character(unlist(snvs$ALT)))
  }
  
  snvs$trinuc <- getSeq(resize(snvs, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
  snvs$ref <- subseq(snvs$trinuc, 2,2)
  flipidx <- which(snvs$ref %in% c("A", "G"))
  snvs$trinuc[flipidx] <- reverseComplement(snvs$trinuc[flipidx])
  snvs$alt <- DNAStringSet(snvs$alt)
  snvs$alt[flipidx] <- reverseComplement(snvs$alt[flipidx])
  
  return(snvs)
}

### remove all driver muts to reduce influence of selection
# debug(get_snvs)
allsnvs <- unlist(GRangesList(mclapply(X = sampleids, FUN = get_snvs, snvsmnvdir = SNVMNVINDELDIR, mc.preschedule = T, mc.cores = 16)))
allsnvs <- c(allsnvs, thirdhitsforaddition)
allsnvs$trinucmut <- factor(paste0(allsnvs$trinuc, ">", allsnvs$alt))
allsnvs <- subsetByOverlaps(x = allsnvs, ranges = driverhits, invert = T)

# saveRDS(object = allsnvs, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/mutability_check_allsnvs.rds")

# allsnvs$trinuc <- getSeq(resize(allsnvs, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
# allsnvs$ref <- subseq(allsnvs$trinuc, 2,2)
# flipidx <- which(allsnvs$ref %in% c("A", "G"))
# allsnvs$trinuc[flipidx] <- reverseComplement(allsnvs$trinuc[flipidx])
# allsnvs$alt <- DNAStringSet(allsnvs$alt)
# allsnvs$alt[flipidx] <- reverseComplement(allsnvs$alt[flipidx])
# allsnvs$trinucmut <- factor(paste0(allsnvs$trinuc, ">", allsnvs$alt))

# do per trinuc type
allsnvs_dedup <- sort(unique(allsnvs))
mcols(allsnvs_dedup) <- DataFrame(trinuc = allsnvs_dedup$trinuc)
allsnvs_dedup$is_parhit <- overlapsAny(query = allsnvs_dedup, subject = parhits, type = "equal") | overlapsAny(query = allsnvs_dedup, subject = thirdhitsforaddition, type = "equal")
allsnvs_dedup$hitcount <- countOverlaps(query = allsnvs_dedup, subject = allsnvs, type = "equal") - 1
allsnvs_dedup <- split(x = allsnvs_dedup, f = as.character(allsnvs_dedup$trinuc))


get_mutation_rate <- function(allsnvs_dedup, trinuc, nsim = 100) {
  print(paste0("Permutation testing ", trinuc))
# trinuc <- 'TCT'
  typed_dedup_muts <- allsnvs_dedup[[trinuc]]
  obs <- aggregate(x = typed_dedup_muts$hitcount, by = list(parhit = typed_dedup_muts$is_parhit), FUN = mean)
  
  estim <- replicate(nsim, expr = aggregate(x = typed_dedup_muts$hitcount,
                                            by = list(parhit = typed_dedup_muts$is_parhit),
                                            FUN = function(x) mean(x = sample(x = x, size = length(x), replace = T))), simplify = F)
  quants <- quantile(x = sapply(X = estim, FUN = function(x) x$x[2]/x$x[1]), probs = c(0.025, 0.5, 0.975))
  
  sim <- replicate(nsim, expr = aggregate(x = typed_dedup_muts$hitcount,
                                        by = list(parhit = sample(typed_dedup_muts$is_parhit,size = length(typed_dedup_muts), replace = F)),
                                        FUN = mean), simplify = F)
  pval <- sum(obs$x[2]/obs$x[1] <= sapply(X = sim, FUN = function(x) x$x[2]/x$x[1]))/nsim
  
  outdf <- data.frame(trinuc = trinuc, rate_single = obs$x[1], rate_par = obs$x[2], pval = pval, npar = sum(typed_dedup_muts$is_parhit), ci_lo = quants[[1]], ci_med = quants[[2]], ci_hi = quants[[3]])
  return(outdf)
}

# outrates <- mclapply(X = allsnvs_dedup, FUN = function(x) by(data = x$hitcount, INDICES = x$is_parhit, table), mc.preschedule = T, mc.cores = 10)

# debug(get_mutation_rate)
outrates <- do.call(rbind, mclapply(X = names(allsnvs_dedup), FUN = get_mutation_rate, allsnvs_dedup = allsnvs_dedup, nsim = 1000, mc.preschedule = T, mc.cores = 16))

write.table(x = outrates, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/mutation_rate_comparison_at_par-nonpar_loci.txt", quote = F, col.names = T, row.names = F, sep = "\t")

# 
### some checks on output

library(ggplot2)
outrates <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/mutation_rate_comparison_at_par-nonpar_loci.txt", as.is = T)
outrates2 <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/mutation_rate_comparison_at_par-nonpar_loci_typed.txt", as.is = T)

outrates$trinuc <- factor(outrates$trinuc, levels = outrates$trinuc[order(outrates$npar, decreasing = T)])
p1 <- ggplot(data = outrates, mapping = aes(x = trinuc)) + geom_linerange(mapping = aes(ymin = ci_lo, ymax = ci_hi)) + geom_point(mapping = aes(y = rate_par/rate_single, size = log10(npar)), colour = "red", alpha = .7)
p1 <- p1 + theme_minimal() + scale_y_log10(breaks = c(1,10,100), limits = c(1, 110)) + annotation_logticks(sides = "l") + theme(axis.text.x = element_text(angle = 90))
p1 <- p1 + labs(y = "Ratio mu_parallel/mu_single", x = "")
p1
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/mutation_rate_comparison_at_par-nonpar_loci.pdf"), plot = p1, width = 8, height = 5, useDingbats=FALSE)


p2 <- ggplot(data = outrates2, mapping = aes(x = trinucmut)) + geom_pointrange(mapping = aes(y = ci_med, ymin = ci_lo, ymax = ci_hi)) + geom_point(mapping = aes(y = rate_par/rate_single, size = log10(npar)), colour = "red", alpha = .7)
p2 <- p2 + scale_y_log10() + theme(axis.text.x = element_text(angle = 90))
p2

library(gridExtra)
grid.arrange(p1, p2, nrow = 2)
