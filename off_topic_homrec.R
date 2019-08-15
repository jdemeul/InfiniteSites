### quick check on samples with lots of cnLOH
library(parallel)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(VariantAnnotation)
library(Exact)

tgl <- sum(seqlengths(BSgenome.Hsapiens.1000genomes.hs37d5)[1:24])

SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

samplecnfiles <- list.files(path = "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/",
                            pattern = "*.consensus.20170119.somatic.cna.annotated.txt", full.names = T)
samplecns <- mclapply(X = samplecnfiles, FUN = read.delim, as.is = T, mc.preschedule = T, mc.cores = 10)
samples <- sub(pattern = ".consensus.20170119.somatic.cna.annotated.txt", fixed = T, replacement = "", x = basename(samplecnfiles))
names(samplecns) <- samples

outstats <- mclapply(X = samplecns, FUN = function(df) {
  diplength <- sum((df$end-df$start)[which(df$total_cn == 2, )], na.rm = T)
  cnlohlength <- sum((df$end-df$start)[which(df$major_cn == 2 & df$minor_cn == 0, )], na.rm = T)
  homlength <- sum((df$end-df$start)[which(df$minor_cn == 0, )], na.rm = T)
  return(c(diplength = diplength, cnlohlength = cnlohlength, homlength = homlength))
}, mc.preschedule = T, mc.cores = 10)

outstats <- data.frame(do.call(rbind, outstats), sampleid = samples)

p1 <- ggplot(outstats) + geom_point(mapping = aes(x = diplength, y = cnlohlength))
p1

p1 <- ggplot(outstats) + geom_histogram(mapping = aes(x = cnlohlength/diplength)) + scale_y_log10()
p1

outstats$histology <- sumtable_whole[match(x = outstats$sampleid, table = sumtable_whole$samplename), "histology_abbreviation"]
outstats$ploidy <- sumtable_whole[match(x = outstats$sampleid, table = sumtable_whole$samplename), "ploidy"]
outstats$frachomo <- outstats$homlength/tgl
outstats$frachomoall <- outstats$cnlohlength/tgl




#### check CN profiles for steps with total CN == & minor CN !=

temp <- samplecns[[1]]

potxocn <- lapply(samplecns, FUN = function(df) {
  potxo <- c(df$chromosome[-nrow(df)] == df$chromosome[-1] &
              df$total_cn[-nrow(df)] == df$total_cn[-1] &
              df$minor_cn[-nrow(df)] != df$minor_cn[-1] &
              df$level[-nrow(df)] %in% letters[1:4] & df$level[-1] %in% letters[1:4])
  df$potxo <- c(c(potxo, F) | c(F, potxo))
  return(df[which(df$potxo),])
})

potxocndf <- do.call(rbind, potxocn)
potxocndf$sampleid <- rep(names(potxocn), sapply(potxocn, nrow))

sumtable_whole$npotxo <- c(table(potxocndf$sampleid))[sumtable_whole$samplename]

p1 <- ggplot(data = sumtable_whole, mapping = aes(x = npotxo, y = ploidy)) + geom_point()
p1



### proving mitotic recombination
# 1) get all steps with total CN == & minor CN !=
# check absence of SV
# 2) get all sites where phasing normal-tumour is inconsitent
# 3) overlap within 10kb? (will be many inconsistent in high mut load tumours)

sampleid <- "a29278af-7ecf-403e-b6a9-623ea7879d05"
CNDIR <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"
PHASDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20181021_hetSNPs_all+phasing_out/"
SVDIR <- "/srv/shared/vanloo/ICGC-structural-variants/"

get_xo_info <- function(sampleid, cndir, phasdir, svdir) {

  print(paste0("Running sample ", sampleid))
  
# read data
cnfile <- file.path(cndir, paste0(sampleid, ".consensus.20170119.somatic.cna.annotated.txt"))
phasfilet <- file.path(phasdir, sampleid, paste0(sampleid, "_tumour_snv-snp_phased.txt"))
phasfilen <- file.path(phasdir, sampleid, paste0(sampleid, "_normal_snv-snp_phased.txt"))
svfile <- file.path(svdir, paste0(sampleid, ".pcawg_consensus_1.6.161116.somatic.sv.vcf.gz"))

if (!all(file.exists(cnfile, phasfilen, phasfilen, svfile))) return(NULL)

cndata <- read.delim(file = cnfile, as.is = T)
potxo <- c(cndata$chromosome[-nrow(cndata)] == cndata$chromosome[-1] &
             cndata$total_cn[-nrow(cndata)] == cndata$total_cn[-1] &
             cndata$minor_cn[-nrow(cndata)] != cndata$minor_cn[-1] &
             cndata$level[-nrow(cndata)] %in% letters[1:4] & cndata$level[-1] %in% letters[1:4])
cndata <- cndata[which(c(c(potxo, F) | c(F, potxo))), ]

if (nrow(cndata) == 0) return(NULL)

relsegs <- reduce(GRanges(cndata, seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
relsegs <- resize(x = relsegs, width = width(relsegs) - 2, fix = "center")
relbps <- subsetByOverlaps(x = GRanges(seqnames = cndata$chromosome, ranges = IRanges(start = cndata$start, end = cndata$start), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)),
                           ranges = relsegs)

svs <- rowRanges(readVcf(file = svfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
# relsvs <- subsetByOverlaps(x = svs, ranges = relbps)
relbps$isconsSV <- overlapsAny(query = relbps, subject = svs)

phasn <- read.delim(file = phasfilen, as.is = T)
phast <- read.delim(file = phasfilet, as.is = T)

phasall <- merge(x = phasn, y = phast, by = c("chr", "pos1", "ref1", "alt1", "type1", "pos2", "ref2", "alt2", "type2"))
colnames(phasall) <- sub(pattern = ".x", replacement = "_norm", fixed = T, x = sub(pattern = ".y", replacement = "_tum", x = colnames(phasall), fixed = T))
phasall <- phasall[which(phasall$type1 == "SNP" & phasall$type2 == "SNP"), ]

phasall <- phasall[overlapsAny(query = GRanges(seqnames = phasall$chr, ranges = IRanges(start = phasall$pos1, end = phasall$pos2), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)),
            subject = resize(relbps, width = ifelse(relbps$isconsSV, 1, 2e5), fix = "center")), ]

testmat <- cbind(phased_norm = rowSums(phasall[, c("Num_ref_ref_norm", "Num_alt_alt_norm")]), 
                 anti_norm = rowSums(phasall[, c("Num_ref_alt_norm", "Num_alt_ref_norm")]),
                 phased_tum = rowSums(phasall[, c("Num_ref_ref_tum", "Num_alt_alt_tum")]),
                 anti_tum = rowSums(phasall[, c("Num_ref_alt_tum", "Num_alt_ref_tum")]))
phasall$p_imbalance <- apply(X = testmat, MARGIN = 1, FUN = function(x) {
  singlemat <- matrix(x, byrow = T, nrow = 2)
  if (any(c(rowSums(singlemat), colSums(singlemat)) == 0)) 
    return(1)
  else
    exact.test(data = singlemat, alternative = "two.sided", to.plot = F, method = "Boschloo")$p.value
  })

return(list(breaks = cndata, bps = relbps, phasing = phasall))
}


SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

# get_xo_info(sampleid = sampleid, cndir = CNDIR, phasdir = PHASDIR, svdir = SVDIR)

outlist <- mclapply(sumtable_whole$samplename, FUN = get_xo_info, cndir = CNDIR, phasdir = PHASDIR, svdir = SVDIR, mc.preschedule = T, mc.cores = 16)
names(outlist) <- sumtable_whole$samplename
saveRDS(object = outlist, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/off_topic_mitXO.RDS")


hasphaserr <- sapply(X = outlist, FUN = function(x) any(x$phasing$p_imbalance < .05))
names(hasphaserr) <- sumtable_whole$samplename

which(hasphaserr)
# 45d0ccb2-641f-4348-b3a8-61f4113cd85b is canadian prostate CA sample which was not flagged up (but has quite a bit of converted stuff)
# 8c233a11-3b2e-4273-bbe1-b5a5f5a351d5 is PACA-CA (also lot of weird converted)

# interesting:
# a7fb0931-28df-46f4-bc0f-2011fc91f0e1 (PBCA)
# b7722577-f200-4dec-97f6-4fab3a9ba52d (CLL)
# ba94c29b-b76e-4d67-bf5a-ce6bc45d85f8 (LGG) CNS-Oligo (also interesting region with switches back and forth)
# ec43c4b5-fb72-4a4a-af03-10c2d05ff159 (COAD)

# unsure
# bb567851-d4ff-4a93-8576-04a37aea68af (RECA)
# d2ab2555-7288-47a4-a80c-bf62d65b67b8 (LUSC)



outlist[["ec43c4b5-fb72-4a4a-af03-10c2d05ff159"]]

View(outlist[["ec43c4b5-fb72-4a4a-af03-10c2d05ff159"]]$phasing)


# check # of breaks with/o SV
sum(sapply(X = outlist, FUN = function(x) sum(x$bps$isconsSV), simplify = T))
# = 283
sum(sapply(X = outlist, FUN = function(x) length(x$bps), simplify = T))
# = 735


get_xo_info_ctrl <- function(sampleid, cndir, phasdir, svdir) {
  
  print(paste0("Running sample ", sampleid))
  
  # read data
  cnfile <- file.path(cndir, paste0(sampleid, ".consensus.20170119.somatic.cna.annotated.txt"))
  phasfilet <- file.path(phasdir, sampleid, paste0(sampleid, "_tumour_snv-snp_phased.txt"))
  phasfilen <- file.path(phasdir, sampleid, paste0(sampleid, "_normal_snv-snp_phased.txt"))
  svfile <- file.path(svdir, paste0(sampleid, ".pcawg_consensus_1.6.161116.somatic.sv.vcf.gz"))
  
  if (!all(file.exists(cnfile, phasfilen, phasfilen, svfile))) return(NULL)
  
  cndata <- read.delim(file = cnfile, as.is = T)
  potxo <- c(cndata$chromosome[-nrow(cndata)] == cndata$chromosome[-1] &
               (cndata$total_cn[-nrow(cndata)] > 1 | cndata$total_cn[-1] > 1) &
               xor(x = cndata$major_cn[-nrow(cndata)] != cndata$major_cn[-1], cndata$minor_cn[-nrow(cndata)] != cndata$minor_cn[-1]) &
               cndata$level[-nrow(cndata)] %in% letters[1:4] & cndata$level[-1] %in% letters[1:4])
  cndata <- cndata[which(c(c(potxo, F) | c(F, potxo))), ]
  
  if (nrow(cndata) == 0) return(NULL)
  
  relsegs <- reduce(GRanges(cndata, seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
  relsegs <- resize(x = relsegs, width = width(relsegs) - 2, fix = "center")
  relbps <- subsetByOverlaps(x = GRanges(seqnames = cndata$chromosome, ranges = IRanges(start = cndata$start, end = cndata$start), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)),
                             ranges = relsegs)
  
  svs <- rowRanges(readVcf(file = svfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
  # relsvs <- subsetByOverlaps(x = svs, ranges = relbps)
  relbps$isconsSV <- overlapsAny(query = relbps, subject = svs)
  
  # phasn <- read.delim(file = phasfilen, as.is = T)
  # phast <- read.delim(file = phasfilet, as.is = T)
  # 
  # phasall <- merge(x = phasn, y = phast, by = c("chr", "pos1", "ref1", "alt1", "type1", "pos2", "ref2", "alt2", "type2"))
  # colnames(phasall) <- sub(pattern = ".x", replacement = "_norm", fixed = T, x = sub(pattern = ".y", replacement = "_tum", x = colnames(phasall), fixed = T))
  # phasall <- phasall[which(phasall$type1 == "SNP" & phasall$type2 == "SNP"), ]
  # 
  # phasall <- phasall[overlapsAny(query = GRanges(seqnames = phasall$chr, ranges = IRanges(start = phasall$pos1, end = phasall$pos2), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)),
  #                                subject = resize(relbps, width = ifelse(relbps$isconsSV, 1, 2e5), fix = "center")), ]
  # 
  # testmat <- cbind(phased_norm = rowSums(phasall[, c("Num_ref_ref_norm", "Num_alt_alt_norm")]), 
  #                  anti_norm = rowSums(phasall[, c("Num_ref_alt_norm", "Num_alt_ref_norm")]),
  #                  phased_tum = rowSums(phasall[, c("Num_ref_ref_tum", "Num_alt_alt_tum")]),
  #                  anti_tum = rowSums(phasall[, c("Num_ref_alt_tum", "Num_alt_ref_tum")]))
  # phasall$p_imbalance <- apply(X = testmat, MARGIN = 1, FUN = function(x) {
  #   singlemat <- matrix(x, byrow = T, nrow = 2)
  #   if (any(c(rowSums(singlemat), colSums(singlemat)) == 0)) 
  #     return(1)
  #   else
  #     exact.test(data = singlemat, alternative = "two.sided", to.plot = F, method = "Boschloo")$p.value
  # })
  
  return(list(breaks = cndata, bps = relbps, phasing = NULL))
}


ctrllist <- mclapply(sumtable_whole$samplename, FUN = get_xo_info_ctrl, cndir = CNDIR, phasdir = PHASDIR, svdir = SVDIR, mc.preschedule = T, mc.cores = 16)
names(ctrllist) <- sumtable_whole$samplename


# check # of breaks with/o SV
sum(sapply(X = outlist, FUN = function(x) sum(x$bps$isconsSV), simplify = T))
# = 283
sum(sapply(X = outlist, FUN = function(x) length(x$bps), simplify = T))
# = 735


# check # of breaks with/o SV
sum(sapply(X = ctrllist, FUN = function(x) sum(x$bps$isconsSV), simplify = T))
# = 148569
sum(sapply(X = ctrllist, FUN = function(x) length(x$bps), simplify = T))
# = 171755


df <- data.frame(sampleid = sumtable_whole$samplename, histology = sumtable_whole$histology_abbreviation, is_pref = sumtable_whole$is_preferred,
                 nxo = sapply(X = outlist, FUN = function(x) length(x$bps), simplify = T),
                 nbp = sapply(X = ctrllist, FUN = function(x) length(x$bps), simplify = T))
df$enrich <- df$nxo / (df$nxo+df$nbp)

df$histology <- factor(df$histology)
sort(table(df[df$nxo > 0, "histology"]) / table(df[df$nxo + df$nbp > 0, "histology"]))
sort(table(df[df$nxo > 0, "histology"]) / table(df[, "histology"]))
sort(c(by(data = df$nxo, INDICES = df$histology, FUN = sum)) / c(by(data = df$nbp+df$nxo, INDICES = df$histology, FUN = sum)))
sort(c(by(data = df$nxo, INDICES = df$histology, FUN = sum)) / c(by(data = df$nbp+df$nxo, INDICES = df$histology, FUN = length)))
saveRDS(object = ctrllist, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/off_topic_mitXO_ctrl.RDS")
