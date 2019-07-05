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
parhits$trinuc <- getSeq(resize(parhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
flipparidx <- which(parhits$ref %in% c("A", "G"))
parhits$trinuc[flipparidx] <- reverseComplement(parhits$trinuc[flipparidx])
parhits$alt <- DNAStringSet(parhits$alt)
parhits$alt[flipparidx] <- reverseComplement(parhits$alt[flipparidx])
parhits$trinucmut <- paste0(parhits$trinuc, ">", parhits$alt)

# thirdhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T))
# thirdhits <- thirdhits[!seqnames(thirdhits) == "Y"]

# thirdhitsforaddition <- granges(thirdhits)
# mcols(thirdhitsforaddition) <- DataFrame(alt = thirdhits$alt1,
#                                          trinuc = getSeq(resize(thirdhitsforaddition, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5),
#                                          ref = thirdhits$REF)
# flipthirdidx <- which(thirdhitsforaddition$ref %in% c("A", "G"))
# thirdhitsforaddition$trinuc[flipthirdidx] <- reverseComplement(thirdhitsforaddition$trinuc[flipthirdidx])
# thirdhitsforaddition$alt <- DNAStringSet(thirdhitsforaddition$alt)
# thirdhitsforaddition$alt[flipthirdidx] <- reverseComplement(thirdhitsforaddition$alt[flipthirdidx])


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
# allsnvs <- c(allsnvs, thirdhitsforaddition)
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


get_mutation_rate <- function(allsnvs, trinucmut, parhits, nsim = 100) {
  print(paste0("Permutation testing ", trinucmut))
  trinuc <- substring(trinucmut,1,3)
  
  typedidxs <- which(allsnvs$trinucmut == trinucmut)
  typed_dedup_muts <- sort(unique(allsnvs[typedidxs]))
  
  typed_dedup_muts$hitcount <- countOverlaps(query = typed_dedup_muts,
                                             subject = allsnvs[typedidxs],
                                             type = "equal") - 1
  # trinuc <- 'TCT'
  typed_dedup_muts$is_parhit <- overlapsAny(query = typed_dedup_muts, subject = parhits[which(parhits$trinucmut == trinucmut)], type = "equal")
  obs <- aggregate(x = typed_dedup_muts$hitcount, by = list(parhit = typed_dedup_muts$is_parhit), FUN = mean)
  
  estim <- replicate(nsim, expr = aggregate(x = typed_dedup_muts$hitcount,
                                            by = list(parhit = typed_dedup_muts$is_parhit),
                                            FUN = function(x) mean(x = sample(x = x, size = length(x), replace = T))), simplify = F)
  quants <- quantile(x = sapply(X = estim, FUN = function(x) x$x[2]/x$x[1]), probs = c(0.025, 0.5, 0.975))
  
  sim <- replicate(nsim, expr = aggregate(x = typed_dedup_muts$hitcount,
                                          by = list(parhit = sample(typed_dedup_muts$is_parhit,size = length(typed_dedup_muts), replace = F)),
                                          FUN = mean), simplify = F)
  pval <- sum(obs$x[2]/obs$x[1] <= sapply(X = sim, FUN = function(x) x$x[2]/x$x[1]))/nsim
  
  outdf <- data.frame(trinucmut = trinucmut, rate_single = obs$x[1], rate_par = obs$x[2], pval = pval, npar = sum(typed_dedup_muts$is_parhit), ci_lo = quants[[1]], ci_med = quants[[2]], ci_hi = quants[[3]])
  return(outdf)
}

# outrates <- mclapply(X = allsnvs_dedup, FUN = function(x) by(data = x$hitcount, INDICES = x$is_parhit, table), mc.preschedule = T, mc.cores = 10)

# debug(get_mutation_rate)
outrates <- do.call(rbind, mclapply(X = levels(allsnvs$trinucmut), FUN = get_mutation_rate, allsnvs = allsnvs, parhits = parhits, nsim = 1000, mc.preschedule = F, mc.cores = 16))

write.table(x = outrates, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/mutation_rate_comparison_at_par-nonpar_loci_typed.txt", quote = F, col.names = T, row.names = F, sep = "\t")

# 
# 
# rateTTT <- outrates$TCT
# quantile(sample(x = as.integer(names(rateTTT$'TRUE')), prob = rateTTT$'TRUE', size = 1e5, replace = T)/
#   sample(x = as.integer(names(rateTTT$'FALSE')), prob = rateTTT$'FALSE', size = 1e5, replace = T), probs = c(0.025, 0.5, 0.975))
# 
# temp <- as.data.frame(mcols(allsnvs_dedup)[which(allsnvs_dedup$trinuc == "TCT"),])
# # temp <- by(data = allsnvs_dedup$hitcount, INDICES = allsnvs_dedup$trinuc, FUN = table)
# 
# 
# allsnvoverlaps <- findOverlaps(query = allsnvs_dedup, subject = allsnvs, type = "equal")
# 
# # temp <- do.call(rbind, (by(data = subjectHits(allsnvoverlaps), INDICES = queryHits(allsnvoverlaps), FUN = function(subidxs, allsnvs) table(allsnvs$alt[subidxs]),
# #    allsnvs = allsnvs)))
# 
# allsnvoverlaps$hittypes
# 
# # snvs$sampleid <- rep(sampleid, length(snvs))
# 
# # annotate with trinuc context
# sampleranges$trinuc <- getSeq(sampleranges, x = BSgenome.Hsapiens.1000genomes.hs37d5)
# agidxs <- which(as.character(subseq(sampleranges$trinuc, 2, 2)) %in% c("A", "G"))
# mcols(sampleranges)[agidxs, "trinuc"] <- reverseComplement(mcols(sampleranges)[agidxs, "trinuc"])
# 
# genometiles <- unlist(tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.1000genomes.hs37d5)[1], tilewidth = 1))
# 
# 
