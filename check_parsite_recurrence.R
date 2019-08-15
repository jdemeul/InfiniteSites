### checking recurrent biallelics


# load all variants
library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
# library(ggplot2)
library(VariantAnnotation)
library(parallel)
# library(rslurm)
# library(VGAM)
library(ggplot2)

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
parhits <- parhits[parhits$sampleid %in% sampleids]
thirdhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T))
thirdhits <- thirdhits[!seqnames(thirdhits) == "Y"]
thirdhits <- thirdhits[thirdhits$sampleid %in% sampleids]

thirdhitsforaddition <- granges(thirdhits)
mcols(thirdhitsforaddition) <- DataFrame(alt = thirdhits$alt1,
                                         trinuc = getSeq(resize(thirdhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5),
                                         ref = thirdhits$REF)
# thirdhitsforaddition <- c(granges(thirdhits),granges(thirdhits))
# mcols(thirdhitsforaddition) <- DataFrame(alt = c(thirdhits$alt1,thirdhits$alt2),
#                                          trinuc = rep(getSeq(resize(thirdhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5), 2),
#                                          ref = rep(thirdhits$REF, 2))
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
  snvs$sampleid <- rep(sampleid, length(snvs))
  # snvs <- snvs[which(snvs$trinuc == "TCT")]
  
  return(snvs)
}

### remove all driver muts to reduce influence of selection
# debug(get_snvs)
allsnvs <- unlist(GRangesList(mclapply(X = sampleids, FUN = get_snvs, snvsmnvdir = SNVMNVINDELDIR, mc.preschedule = T, mc.cores = 10)))
# allsnvs <- c(allsnvs, thirdhitsforaddition)#[which(thirdhitsforaddition$trinuc == "TCT")])
# allsnvs$trinucmut <- factor(paste0(allsnvs$trinuc, ">", allsnvs$alt))
# allsnvs <- subsetByOverlaps(x = allsnvs, ranges = driverhits, invert = T)
# 
# # saveRDS(object = allsnvs, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/mutability_check_allsnvs+bial.rds")
# allsnvs <- readRDS(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/mutability_check_allsnvs+bial.rds")
# 
# 
# allsnvs_dedup <- granges(x = sort(unique(allsnvs)), use.mcols = T, use.names = F)
# # mcols(allsnvs_dedup) <- DataFrame(trinuc = allsnvs_dedup$trinuc)
# # colnames(mcols(allsnvs_dedup)) <- "trinuc"
# 
# # allsnvs <- granges(allsnvs, use.names=FALSE, use.mcols=FALSE)
# 
# allsnvs_dedup$is_parhit <- as.factor(overlapsAny(query = allsnvs_dedup, subject = c(granges(parhits), granges(thirdhits)), type = "equal"))
# allsnvs_dedup$hitcount <- countOverlaps(query = allsnvs_dedup, subject = allsnvs, type = "equal")
# allsnvs_dedup$trinuc <- as.factor(allsnvs_dedup$trinuc)
# allsnvs_dedup_dfonly <- mcols(allsnvs_dedup, use.names = F)[, c("trinuc", "is_parhit", "hitcount")]
# 
# 
# 


#### get all loci with par or bial variants

bialvar <- c(granges(parhits, use.names = F, use.mcols = F), granges(thirdhits, use.names = F, use.mcols = F))
bialloci <- unique(bialvar)
mcols(bialloci)$nbial <- countOverlaps(query = bialloci, subject = bialvar, type = "equal")
table(bialloci$nbial)

bialloci[bialloci$nbial == 4]
hitpars <- subsetByOverlaps(x = parhits, ranges = bialloci[bialloci$nbial == 4], type = "equal")
hitthird <- subsetByOverlaps(x = thirdhits, ranges = bialloci[bialloci$nbial == 4], type = "equal")
hitsnvs <- subsetByOverlaps(x = allsnvs, ranges = bialloci[bialloci$nbial == 4], type = "equal")

hitsnvsdf <- as.data.frame(mcols(hitsnvs))
hitsnvsdf$type <- ifelse(hitsnvsdf$sampleid %in% hitpars$sampleid, "par", "mono")
hitsnvsdf <- hitsnvsdf[!hitsnvsdf$sampleid %in% hitthird$sampleid, ]
hitsnvsdf <- rbind(hitsnvsdf, hitsnvsdf[hitsnvsdf$type == "par", ], data.frame(alt = c(hitthird$alt1, hitthird$alt2), trinuc = rep("TCT", 2), ref = rep("C", 2), sampleid = rep(hitthird$sampleid, 2), type = rep("bial",2)))
hitsnvsdf <- hitsnvsdf[order(hitsnvsdf$sampleid), ]
# top site has
hitsnvsdf$histology_abbreviation <- sumtab[match(x = hitsnvsdf$sampleid, table = sumtab$samplename), "histology_abbreviation"]
hitsnvsdf$chr <- "19"
hitsnvsdf$pos <- 17970682
rownames(hitsnvsdf) <- NULL

write.table(x = hitsnvsdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/TopBialHitSite.txt", quote = F, row.names = F, col.names = T, sep = "\t")

### check how many biallelic SNVs are monoallelically called
calledbials <- subsetByOverlaps(x = allsnvs, ranges = thirdhits, type = "equal")
calledbials <- calledbials[paste0(calledbials$sampleid, "_", start(calledbials)) %in% paste0(thirdhits$sampleid, "_", start(thirdhits))]




# check DBS parallel variants
bialloci[overlapsAny(query = bialloci, drop.self = T, maxgap = 0)]

dbstretches <- reduce(x = parhits)
dbstretches <- dbstretches[width(dbstretches) == 2]
pardbs <- subsetByOverlaps(x = parhits, ranges = dbstretches)

temp <- unlist(GRangesList(lapply(X = split(parhits, parhits$sampleid), FUN = function(x) x[overlapsAny(query = x, drop.self = T, maxgap = 0)])))
temp$heptanuc <- getSeq(resize(temp, fix = "center", width = 15), x = BSgenome.Hsapiens.1000genomes.hs37d5)
names(temp) <- NULL
View(as.data.frame(temp))

write.table(x = as.data.frame(temp), file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/Parallel_dinucleotidevariants.txt", quote = F, row.names = F, col.names = T, sep = "\t")

