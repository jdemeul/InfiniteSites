### checking new simulation results
# eg sample 93ff... with ~ 1M variants

library(VariantAnnotation)
library(parallel)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

# snvs <- readVcf(file = "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/93ff786e-0165-4b02-8d27-806d422e93fc.consensus.20160830.somatic.snv_mnv.vcf.gz", genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
snvs <- readVcf(file = "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/14c5b81d-da49-4db1-9834-77711c2b1d38.consensus.20160830.somatic.snv_mnv.vcf.gz", genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))


multihits <- findOverlaps(snvs, drop.self = T, drop.redundant = T)
multihits

any(lengths(mcols(snvs)$ALT) > 1)


### check all samples for multi-allele sites

vcffiles <- list.files(path = "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/", pattern = ".consensus.20160830.somatic.snv_mnv.vcf.gz$", full.names = T, recursive = T)
sampleids <- gsub(pattern = ".consensus.20160830.somatic.snv_mnv.vcf.gz", replacement = "", x = basename(vcffiles))

hasmultialleles <- mclapply(X = vcffiles, FUN = function(x) {y <- readVcf(x, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)); y[overlapsAny(y, drop.self = T, drop.redundant = F)]}, mc.cores = 16, mc.preschedule = F)
hasmultialleleslength <- sapply(hasmultialleles, length)

sampleids[which(hasmultialleleslength > 0)]
lapply(hasmultialleles[which(hasmultialleleslength > 0)], info)

any(unlist(hasmultialleles))

lengths(mcols(rowRanges(snvs))$ALT)

mutectsnvs <- readVcf("/srv/shared/vanloo/ICGC/READ-US/SNVs_preconcensus/broad_variant_calling_pipeline/84c07c54-62c6-4a3b-811c-551adc1f5d52/14c5b81d-da49-4db1-9834-77711c2b1d38.broad-mutect-v3.20160222.somatic.snv_mnv.vcf.gz", genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
templengths <- lengths(rowRanges(mutectsnvs)$ALT)
any(templengths > 1)
temp <- rowRanges(mutectsnvs)$ALT


rowRanges(mutectsnvs)[overlapsAny(mutectsnvs, drop.self = T, drop.redundant = F)]

rowRanges(subsetByOverlaps(x = mutectsnvs, ranges = GRanges(seqnames = "10", ranges = IRanges(start = 68727788, end = 68727790))))

mutectsnvs <- readVcf("/srv/shared/vanloo/ICGC/READ-US/SNVs_preconcensus/dkfz_embl_variant_calling_pipeline/d469e725-fefe-4515-9b65-42584a314201/14c5b81d-da49-4db1-9834-77711c2b1d38.dkfz-snvCalling_1-0-132-1-hpc.1508061844.germline.snv_mnv.vcf.gz", genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
rowRanges(mutectsnvs)[overlapsAny(mutectsnvs, drop.self = T, drop.redundant = F)]
rowRanges(mutectsnvs)[lengths(rowRanges(mutectsnvs)$ALT) > 1]

# proper look at biallelic variants
# sampleid <- "14c5b81d-da49-4db1-9834-77711c2b1d38"
sampleid <- "93ff786e-0165-4b02-8d27-806d422e93fc"
CONSFILE <- paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
preconsensusfiles <- list.files(path = "/srv/shared/vanloo/ICGC/", pattern = paste0("^", sampleid, ".*.somatic.snv_mnv.vcf.gz$"), full.names = T, recursive = T)
BROADFILE <- grep(pattern = "broad", x = preconsensusfiles, value = T)
MUSEFILE <- grep(pattern = "muse", x = preconsensusfiles, value = T)
DKFZFILE <- grep(pattern = "dkfz", x = preconsensusfiles, value = T)
SANGERFILE <- grep(pattern = "sanger", x = preconsensusfiles, value = T)
# BROADFILE <- "/srv/shared/vanloo/ICGC/READ-US/SNVs_preconcensus/broad_variant_calling_pipeline/84c07c54-62c6-4a3b-811c-551adc1f5d52/14c5b81d-da49-4db1-9834-77711c2b1d38.broad-mutect-v3.20160222.somatic.snv_mnv.vcf.gz"
# MUSEFILE <- "/srv/shared/vanloo/ICGC/READ-US/SNVs_preconcensus/muse_variant_calling_pipeline/cde5254c-0e9d-4e8a-ad68-ce8bb13b12a3/14c5b81d-da49-4db1-9834-77711c2b1d38.MUSE_1-0rc-vcf.20150918.somatic.snv_mnv.vcf.gz"
# DKFZFILE <- "/srv/shared/vanloo/ICGC/READ-US/SNVs_preconcensus/dkfz_embl_variant_calling_pipeline/d469e725-fefe-4515-9b65-42584a314201/14c5b81d-da49-4db1-9834-77711c2b1d38.dkfz-snvCalling_1-0-132-1-hpc.1508061844.somatic.snv_mnv.vcf.gz"
# SANGERFILE <- "/srv/shared/vanloo/ICGC/READ-US/SNVs_preconcensus/sanger_variant_calling_pipeline/c6352790-cff9-11e5-96ab-eaf5899203f5/14c5b81d-da49-4db1-9834-77711c2b1d38.svcp_1-0-5.20150403.somatic.snv_mnv.vcf.gz"
consensus <- readVcf(file = CONSFILE, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
broad <- readVcf(file = BROADFILE, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
muse <- readVcf(file = MUSEFILE, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
dkfz <- readVcf(file = DKFZFILE, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
sanger <- readVcf(file = SANGERFILE)

#### Muse only one to regularly call multiallelics, yet these are all hetSNPs where one allele is somatic
#### DKFZ calls one or two, but that's it

# no overlaps
# broadoverlaps <- broad[overlapsAny(query = broad, drop.self = T, drop.redundant = F)]
# broadmal <- broad[lengths(alt(broad)) > 1]
# museoverlaps <- muse[overlapsAny(query = muse, drop.self = T, drop.redundant = F)]
# 1094 calls total
musemal <- muse[lengths(alt(muse)) > 1]
# dkfzoverlaps <- dkfz[overlapsAny(query = dkfz, drop.self = T, drop.redundant = F)]
# 14 calls total
dkfzmal <- dkfz[lengths(alt(dkfz)) > 1]
# sangeroverlaps <- sanger[overlapsAny(query = sanger, drop.self = T, drop.redundant = F)]
# sangermal <- sanger[lengths(alt(sanger)) > 1]

# overlaps
# single snv with two alleles reported as two separate lines
consensusoverlaps <- consensus[overlapsAny(query = consensus, drop.self = T, drop.redundant = F)]
# consensusmal <- consensus[lengths(alt(consensus)) > 1]

deb9ffile <- list.files(path = "/srv/shared/vanloo/ICGC/", pattern = paste0("^93ff786e-0165-4b02-8d27-806d422e93fc.MUSE_1-0rc-vcf.*somatic.snv_mnv.vcf.gz$"), recursive = T, full.names = T)
# muse only one which consistently calls decent number
muse <- readVcf(file = "/srv/shared/vanloo/ICGC//COAD-US/SNVs_preconcensus/dkfz_embl_variant_calling_pipeline/f7998153-1b12-411f-8fd4-759bc8bb8aa1/93ff786e-0165-4b02-8d27-806d422e93fc.dkfz-snvCalling_1-0-132-1-hpc.1508061544.somatic.snv_mnv.vcf.gz", genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))

museclean <- rowRanges(muse[lengths(alt(muse)) > 1 & mcols(muse)$FILTER == "PASS" & seqnames(muse) %in% c(1:22, "X")])
# museclean <- rowRanges(dkfz[lengths(alt(dkfz)) > 1 & seqnames(dkfz) %in% c(1:22, "X")])
mcols(museclean)$trinuc <- getSeq(x = BSgenome.Hsapiens.1000genomes.hs37d5, names = resize(x = museclean, fix = "center", width = 3))
agidxs <- which(as.character(museclean$REF) %in% c("A", "G"))
mcols(museclean)[agidxs, "trinuc"] <- reverseComplement(mcols(museclean)[agidxs, "trinuc"])
mcols(museclean)[agidxs, "REF"] <- reverseComplement(mcols(museclean)[agidxs, "REF"])
mcols(museclean)[agidxs, "ALT"] <- DNAStringSetList(lapply(X = mcols(museclean)[agidxs, "ALT"], FUN = reverseComplement))

# assume all of these have arisen in parallel, i.e. they target distinct REF alleles (neglegible forward mutation)

generate_bases_types_trinuc <- function() {
  # generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  types <- c("A", "G", "T", "A", "C", "G")
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(trinucleotides, ">", rep(types, rep(16,6))))
  return(list(bases = bases, types = types, trinucleotides = trinucleotides,
              trinucleotides_mutations = trinucleotides_mutations))
}

trinucs <- generate_bases_types_trinuc()
idxv <- setNames(object = 1:length(trinucs$trinucleotides_mutations), nm = trinucs$trinucleotides_mutations)

# musemutdf <- data.frame(allele1 = paste0(mcols(museclean)$trinuc, ">", sapply(X = mcols(museclean)$ALT, FUN = function(x) as.character(x[[1]]))), 
#                         allele2 = paste0(mcols(museclean)$trinuc, ">", sapply(X = mcols(museclean)$ALT, FUN = function(x) as.character(x[[2]]))))
mutidx1 <- idxv[paste0(mcols(museclean)$trinuc, ">", sapply(X = mcols(museclean)$ALT, FUN = function(x) as.character(x[[1]])))]
mutidx2 <- idxv[paste0(mcols(museclean)$trinuc, ">", sapply(X = mcols(museclean)$ALT, FUN = function(x) as.character(x[[2]])))]
alhitsdf <- data.frame(mut1 = trinucs$trinucleotides_mutations[pmin(mutidx1, mutidx2)], mut2 = trinucs$trinucleotides_mutations[pmax(mutidx1, mutidx2)])
alhitsdf <- data.frame(table(mut1 = factor(alhitsdf$mut1, levels = trinucs$trinucleotides_mutations), mut2 = factor(alhitsdf$mut2, levels = trinucs$trinucleotides_mutations)))

# annotation
alhitsdf[, c("reftri1", "alt1")] <- do.call(rbind, strsplit(as.character(alhitsdf$mut1), split = ">"))
alhitsdf[, c("reftri2", "alt2")] <- do.call(rbind, strsplit(as.character(alhitsdf$mut2), split = ">"))
alhitsdf$alttri1 <- DNAStringSet(paste0(substr(alhitsdf$reftri1, 1,1), alhitsdf$alt1, substr(alhitsdf$reftri1, 3,3)))
alhitsdf[alhitsdf$alt1 %in% c("A", "G"), "alttri1"] <- reverseComplement(alhitsdf[alhitsdf$alt1 %in% c("A", "G"), "alttri1"])
alhitsdf$alttri1 <- as.character(alhitsdf$alttri1)

alhitsdf$alttri2 <- DNAStringSet(paste0(substr(alhitsdf$reftri2, 1,1), alhitsdf$alt2, substr(alhitsdf$reftri2, 3,3)))
alhitsdf[alhitsdf$alt2 %in% c("A", "G"), "alttri2"] <- reverseComplement(alhitsdf[alhitsdf$alt2 %in% c("A", "G"), "alttri2"])
alhitsdf$alttri2 <- as.character(alhitsdf$alttri2)

alhitsdf$type <- "empty"
alhitsdf$type <- ifelse(alhitsdf$reftri1 == alhitsdf$reftri2,
                      ifelse(alhitsdf$alt1 == alhitsdf$alt2, "parallel", "third_allele"),
                      ifelse(alhitsdf$alttri1 == alhitsdf$reftri2, ifelse(alhitsdf$alttri2 == alhitsdf$reftri1, "back", "forward"), alhitsdf$type))
alhitsdf$type <- factor(alhitsdf$type, levels = c("back", "forward", "parallel", "third_allele", "empty"))
# alhitsdf <- alhitsdf[!(alhitsdf$type == "third_allele" & as.integer(alhitsdf$mut1) > as.integer(alhitsdf$mut2)), ]

library(ggplot2)
p1 <- ggplot(data = alhitsdf[alhitsdf$Freq > 0,], mapping = aes(x = mut1, y = mut2)) + geom_tile(mapping = aes(fill = Freq), show.legend = F)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(drop = F) + scale_y_discrete(drop = F) #+ scale_fill_brewer(type = "qual", palette = "Set1") #+ scale_alpha_discrete(values = c('TRUE' = 1, 'FALSE' = .1))
p1




##### check Mutect2 recalled variants for multiallelics (otherwise PASS variants)
snvs <- readVcf(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/93ff786e-0165-4b02-8d27-806d422e93fc/93ff786e-0165-4b02-8d27-806d422e93fc_mutect2_snvs_indels_oncefiltered.vcf.gz", genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
snvs_germline <- snvs[grepl(pattern = "germline_risk", x = rowRanges(snvs)$FILTER) & lengths(alt(snvs)) == 1]
# identify Mutect2 Germline filter (appears to be P_GERMLINE < -1)
# hist(unlist(info(snvs_germline)$P_GERMLINE))
# get the passed ones
snvs_pass <- snvs[rowRanges(snvs)$FILTER == "PASS"]
# filter for biallelic variants
snvs_pass_biallelic <- snvs_pass[which(lengths(alt(snvs_pass)) == 2)]
# only SNVs
snvs_pass_biallelic <- snvs_pass_biallelic[which(sapply(X = alt(snvs_pass_biallelic), FUN = function(x) all(x %in% c("A", "C", "G", "T"))))]
# both alleles need to have P_GERMLINE < -1, otherwise likely mutation at het SNP
snvs_pass_biallelic <- snvs_pass_biallelic[which(sapply(X = info(snvs_pass_biallelic)$P_GERMLINE, FUN = function(x) all(x < -1)))]

# get read depths
# normal
geno(snvs_pass_biallelic[2])$F1R2[[1]] + geno(snvs_pass_biallelic[2])$F1R2[[1]]
# tumour
geno(snvs_pass_biallelic[2])$F1R2[[2]] + geno(snvs_pass_biallelic[2])$F1R2[[2]]


snvs_clean <- rowRanges(snvs_pass_biallelic)
mcols(snvs_clean)$trinuc <- getSeq(x = BSgenome.Hsapiens.1000genomes.hs37d5, names = resize(x = snvs_clean, fix = "center", width = 3))
agidxs <- which(as.character(snvs_clean$REF) %in% c("A", "G"))
mcols(snvs_clean)[agidxs, "trinuc"] <- reverseComplement(mcols(snvs_clean)[agidxs, "trinuc"])
mcols(snvs_clean)[agidxs, "REF"] <- reverseComplement(mcols(snvs_clean)[agidxs, "REF"])
mcols(snvs_clean)[agidxs, "ALT"] <- DNAStringSetList(lapply(X = mcols(snvs_clean)[agidxs, "ALT"], FUN = reverseComplement))

# assume all of these have arisen in parallel, i.e. they target distinct REF alleles (neglegible forward mutation)

generate_bases_types_trinuc <- function() {
  # generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  types <- c("A", "G", "T", "A", "C", "G")
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(trinucleotides, ">", rep(types, rep(16,6))))
  return(list(bases = bases, types = types, trinucleotides = trinucleotides,
              trinucleotides_mutations = trinucleotides_mutations))
}

trinucs <- generate_bases_types_trinuc()
idxv <- setNames(object = 1:length(trinucs$trinucleotides_mutations), nm = trinucs$trinucleotides_mutations)

# musemutdf <- data.frame(allele1 = paste0(mcols(snvs_clean)$trinuc, ">", sapply(X = mcols(snvs_clean)$ALT, FUN = function(x) as.character(x[[1]]))), 
#                         allele2 = paste0(mcols(snvs_clean)$trinuc, ">", sapply(X = mcols(snvs_clean)$ALT, FUN = function(x) as.character(x[[2]]))))
mutidx1 <- idxv[paste0(mcols(snvs_clean)$trinuc, ">", sapply(X = mcols(snvs_clean)$ALT, FUN = function(x) as.character(x[[1]])))]
mutidx2 <- idxv[paste0(mcols(snvs_clean)$trinuc, ">", sapply(X = mcols(snvs_clean)$ALT, FUN = function(x) as.character(x[[2]])))]
alhitsdf <- data.frame(mut1 = trinucs$trinucleotides_mutations[pmin(mutidx1, mutidx2)], mut2 = trinucs$trinucleotides_mutations[pmax(mutidx1, mutidx2)])
alhitsdf <- data.frame(table(mut1 = factor(alhitsdf$mut1, levels = trinucs$trinucleotides_mutations), mut2 = factor(alhitsdf$mut2, levels = trinucs$trinucleotides_mutations)))

# annotation
alhitsdf[, c("reftri1", "alt1")] <- do.call(rbind, strsplit(as.character(alhitsdf$mut1), split = ">"))
alhitsdf[, c("reftri2", "alt2")] <- do.call(rbind, strsplit(as.character(alhitsdf$mut2), split = ">"))
alhitsdf$alttri1 <- DNAStringSet(paste0(substr(alhitsdf$reftri1, 1,1), alhitsdf$alt1, substr(alhitsdf$reftri1, 3,3)))
alhitsdf[alhitsdf$alt1 %in% c("A", "G"), "alttri1"] <- reverseComplement(alhitsdf[alhitsdf$alt1 %in% c("A", "G"), "alttri1"])
alhitsdf$alttri1 <- as.character(alhitsdf$alttri1)

alhitsdf$alttri2 <- DNAStringSet(paste0(substr(alhitsdf$reftri2, 1,1), alhitsdf$alt2, substr(alhitsdf$reftri2, 3,3)))
alhitsdf[alhitsdf$alt2 %in% c("A", "G"), "alttri2"] <- reverseComplement(alhitsdf[alhitsdf$alt2 %in% c("A", "G"), "alttri2"])
alhitsdf$alttri2 <- as.character(alhitsdf$alttri2)

alhitsdf$type <- "empty"
alhitsdf$type <- ifelse(alhitsdf$reftri1 == alhitsdf$reftri2,
                        ifelse(alhitsdf$alt1 == alhitsdf$alt2, "parallel", "third_allele"),
                        ifelse(alhitsdf$alttri1 == alhitsdf$reftri2, ifelse(alhitsdf$alttri2 == alhitsdf$reftri1, "back", "forward"), alhitsdf$type))
alhitsdf$type <- factor(alhitsdf$type, levels = c("back", "forward", "parallel", "third_allele", "empty"))
# alhitsdf <- alhitsdf[!(alhitsdf$type == "third_allele" & as.integer(alhitsdf$mut1) > as.integer(alhitsdf$mut2)), ]

### check observed vs simulated plot
alhitssim <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown/93ff786e-0165-4b02-8d27-806d422e93fc_infsites_backfwd_allelic.txt", as.is = T)

plotdf <- merge(x = alhitsdf, y = alhitssim, all = F)

# check for forward mutation post-gain gain ~ third allele
plotdf$mut2 <- as.character(plotdf$mut2)
# plotdf$mut2 <- ifelse(plotdf$type == "forward", paste0(plotdf$reftri1, ">", substr(plotdf$alttri2, 2, 2)), plotdf$mut2)
plotdf$mut2 <- ifelse(plotdf$type != "forward", plotdf$mut2,
                                       ifelse(plotdf$alt1 %in% c("A","G"),
                                              paste0(plotdf$reftri1, ">", as.character(reverseComplement(DNAStringSet(plotdf$alt2)))),
                                              paste0(plotdf$reftri1, ">", plotdf$alt2)))

idxv <- setNames(object = 1:length(trinucs$trinucleotides_mutations), nm = trinucs$trinucleotides_mutations)

midx1 <- pmin(idxv[plotdf[plotdf$type == "forward", "mut1"]], idxv[plotdf[plotdf$type == "forward", "mut2"]])
midx2 <- pmax(idxv[plotdf[plotdf$type == "forward", "mut1"]], idxv[plotdf[plotdf$type == "forward", "mut2"]])
plotdf[plotdf$type == "forward", "mut1"] <- trinucs$trinucleotides_mutations[midx1]
plotdf[plotdf$type == "forward", "mut2"] <- trinucs$trinucleotides_mutations[midx2]

plotdf$mut2 <- factor(plotdf$mut2, levels = trinucs$trinucleotides_mutations)

# library(dplyr)
# plotdf %>% filter(type == "forward") %>% group_by(mut1, mut2) %>% summarise(mut1 = )
temp <- split(x = plotdf[plotdf$type == "forward", ], f = paste0(plotdf$mut1, "-", plotdf$mut2)[plotdf$type == "forward"])
temp <- do.call(rbind, lapply(X = temp, FUN = function(x) cbind(x[1, 1:10, drop = F], t(colSums(x[, 11:16, drop = F])))))

plotdf2 <- merge(x = plotdf[plotdf$type == "third_allele", grep(pattern = "freqbf", x = colnames(plotdf), invert = T) ], y = temp[temp$type == "forward", c("mut1", "mut2", "freqbf_low", "freqbf_med", "freqbf_hi", "type")], by = c("mut1", "mut2"))


library(ggplot2)
library(ggrepel)
# p1 <- ggplot(data = plotdf[plotdf$type == "third_allele" & (plotdf$Freq > 0 | plotdf$freqal_hi > 0), ], mapping = aes(x = order(Freq), y = order(freqal_hi))) + geom_jitter()
# p1

# cor(x = plotdf[plotdf$type %in% c("third_allele"), "Freq"], y = plotdf[plotdf$type %in% c("third_allele"), "freqal_hi"], method = "pearson")
# cor(x = plotdf[plotdf$type %in% c("third_allele"), "Freq"], y = plotdf[plotdf$type %in% c("forward"), "freqbf_hi"], method = "pearson")
cor(x = plotdf2$Freq, y = rowSums(plotdf2[, c("freqal_hi", "freqbf_hi")]), method = "pearson")
cor(x = plotdf2$Freq, y = rowSums(plotdf2[, c("freqal_hi", "freqbf_hi")]), method = "spearman")


p1 <- ggplot(data = plotdf2, mapping = aes(x = Freq, y = freqal_hi + freqbf_hi)) + stat_smooth(method = "lm", formula = y ~ 0 + x) 
p1 <- p1 + geom_point(mapping = aes(size = freqal_hi/(freqal_hi+freqbf_hi)), color = "red", alpha = .5, show.legend = T, stroke = 0)
p1 <- p1 + geom_point(mapping = aes(size = freqbf_hi/(freqal_hi+freqbf_hi)), color = "blue", alpha = .5, show.legend = F, stroke = 0)
p1 <- p1 + geom_abline(slope = 1, intercept = 0, alpha =.3) + scale_x_log10() + scale_y_log10() + geom_label_repel(data = plotdf2[plotdf2$freqal_hi + plotdf2$freqbf_hi > 3 | plotdf2$Freq > 3, ], mapping = aes(label = paste0(mut1, " / ", mut2)), alpha = .5)
p1 <- p1 + theme_minimal() + annotation_logticks() + labs(x = "Observed", y = "Simulated") + coord_fixed(ratio = 1, xlim = c(1,1e3), ylim = c(1,1e3))
p1 <- p1 + guides(size = guide_legend(title = "third allele / (third allele + forward)")) + scale_size_area()
p1


