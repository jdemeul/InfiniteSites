### check which hits fall at [C/CC]TTCCG motif

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
library(Biostrings)

# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")
# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")

SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/"
# RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
HISTORELMERGED <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/merged_histology_v6_release_v1.4.txt"

# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
# sampleid <- "fc876f5c-8339-bc9c-e040-11ac0d485160"
# sampleid <- "7cae6c0b-36fe-411b-bbba-093a4c846d84"
# sampleid <- "0332b017-17d5-4083-8fc4-9d6f8fdbbbde"
# sampleid <- "804ffa2e-158b-447d-945c-707684134c87"

sumtab <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
histomerge <- read.delim(file = HISTORELMERGED, as.is = T)
sampleids <- sumtab[sumtab$is_preferred, "samplename"]
sampleids <- intersect(x = sampleids, y = histomerge[histomerge$donor_wgs_included_excluded == "Included", "tumor_wgs_aliquot_id"])

# OLD
# parhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt", as.is = T))
# thirdhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T))
# NEW
parhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_2021.txt", as.is = T))
thirdhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_20210217.txt", as.is = T))


thirdhits$alt <- thirdhits$alt1
thirdhits$ref <- thirdhits$REF
allhits <- c(parhits[, c("ref", "alt", "sampleid")], thirdhits[, c("ref", "alt", "sampleid")])
allhits <- allhits[allhits$sampleid %in% sampleids]
# allhits <- unique(allhits)

allhits$trinuc <- getSeq(resize(allhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
allhits <- allhits[setdiff(x = 1:length(allhits), y = GenomicRanges:::get_out_of_bound_index(resize(allhits, fix = "center", width = 15)))]
allhits$heptanuc <- getSeq(resize(allhits, fix = "center", width = 15), x = BSgenome.Hsapiens.1000genomes.hs37d5)
flipparidx <- which(allhits$ref %in% c("A", "G")) #9
allhits$trinuc[flipparidx] <- reverseComplement(allhits$trinuc[flipparidx])
allhits$heptanuc[flipparidx] <- reverseComplement(allhits$heptanuc[flipparidx])
allhits$ref[flipparidx] <- reverseComplement(DNAStringSet(allhits$ref[flipparidx]))
allhits$alt[flipparidx] <- reverseComplement(DNAStringSet(allhits$alt[flipparidx]))
allhits$histology_abbreviation <- sumtab$histology_abbreviation[match(x = allhits$sampleid, table = sumtab$samplename)]

for (ttype in c("Skin-Melanoma", "ColoRect-AdenoCA", "Eso-AdenoCA")) {
  # ttype <- "Skin-Melanoma"
if (ttype == "Eso-AdenoCA") {
  melahits <- allhits[allhits$histology_abbreviation %in% c(ttype, "Stomach-AdenoCA"), ]
} else {
  melahits <- allhits[allhits$histology_abbreviation == ttype, ]
}
  
# subset to samples with >= 10 biallelic events
  casestolookat <- names(which(table(melahits$sampleid) > 10))
  melahits <- melahits[melahits$sampleid %in% casestolookat, ]
melahits$mut <- paste0(melahits$trinuc, ">", melahits$alt)
# melahits <- allhits[allhits$histology_abbreviation == "ColoRect-AdenoCA", ]
# melahits <- allhits[allhits$histology_abbreviation %in% c("Eso-AdenoCA", "Stomach-AdenoCA"), ]

# View(as.data.frame(melahits)) # look for NFATC4 [TC]TTTC[CA][TA]
# sum(grepl(pattern = "[TC]TTTC[CA][TA]", x = as.character(melahits$heptanuc)))
# sum(grepl(pattern = "[CT]TTCC", x = as.character(melahits$heptanuc)))


# get the set of unmatched sequences: same samples, read all variants and subsample loci which have same muts
melahitids <- unique(melahits$sampleid)

get_snvs <- function(sampleid, snvsmnvdir, bialsnvs, nsamples) {
  print(paste0("Loading snvs for sample ", sampleid))
  snv_mnvfile <- file.path(snvsmnvdir, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz"))
  
  if (!file.exists(snv_mnvfile)) {
    snvs <- GRanges()
    snvs$REF <- DNAStringSet()
    snvs$ALT <- DNAStringSet()
    snvs$sampleid <- sampleid
  } else {
    snvs <- readVcf(file = snv_mnvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
    snvs <- granges(rowRanges(snvs[lengths(alt(snvs)) == 1 & lengths(ref(snvs)) == 1 & seqnames(snvs) != "Y"]), use.names = F, use.mcols = T)
    snvs$ALT <- unlist(snvs$ALT)
    snvs$sampleid <- rep(sampleid, length(snvs))
    # mcols(snvs) <- DataFrame(alt = as.character(unlist(snvs$ALT)))
  }
  
  snvs$trinuc <- getSeq(resize(snvs, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
  if (length(snvs) > 0) 
    snvs <- snvs[setdiff(x = 1:length(snvs), y = GenomicRanges:::get_out_of_bound_index(resize(snvs, fix = "center", width = 15)))]
  snvs$heptanuc <- getSeq(resize(snvs, fix = "center", width = 15), x = BSgenome.Hsapiens.1000genomes.hs37d5)
  snvs <- snvs[grep(pattern = "N", x = snvs$heptanuc, invert = T), ]
  # snvs$ref <- subseq(snvs$trinuc, 2,2)
  flipidx <- which(snvs$REF %in% c("A", "G"))
  snvs$trinuc[flipidx] <- reverseComplement(snvs$trinuc[flipidx])
  snvs$heptanuc[flipidx] <- reverseComplement(snvs$heptanuc[flipidx])
  snvs$REF[flipidx] <- reverseComplement(snvs$REF[flipidx])
  snvs$ALT[flipidx] <- reverseComplement(snvs$ALT[flipidx])
  # snvs <- snvs[which(snvs$trinuc == "TCT")]
  
  # do sampling here
  tosamplemuts <- c(table(mcols(bialsnvs)[bialsnvs$sampleid == sampleid, "mut"]))*nsamples
  snvs$mut <- paste0(snvs$trinuc, ">", snvs$ALT)
  snvs <- snvs[(!overlapsAny(query = snvs, subject = bialsnvs, type = "equal")) & snvs$mut %in% names(tosamplemuts)]
  snvsampled <- unlist(GRangesList(lapply(X = names(tosamplemuts), FUN = function(x, tosample, allsnvs) {
    allsnvs[sample(x = which(allsnvs$mut == x), size = tosample[x], replace = T)]
  }, allsnvs = snvs, tosample = tosamplemuts)))
  
  return(snvsampled)
}

# debug(get_snvs)
# samplesnvs <- get_snvs(sampleid = melahitids[6], snvsmnvdir = SNVMNVINDELDIR, bialsnvs = melahits, nsamples = 10)

# samplesnvs <- lapply(X = melahitids, FUN = get_snvs, snvsmnvdir = SNVMNVINDELDIR)
samplesnvs <- do.call(c, mclapply(X = melahitids, FUN = get_snvs, snvsmnvdir = SNVMNVINDELDIR, bialsnvs = melahits, nsamples = 10, mc.preschedule = T, mc.cores = 5))
# samplesnvs <- samplesnvs[!overlapsAny(query = samplesnvs, drop.self = T)]
# samplesnvs <- subsetByOverlaps(x = samplesnvs, ranges = melahits, type = "equal", invert = T)
# samplesnvs$mut <- paste0(samplesnvs$trinuc, ">", samplesnvs$ALT)



#### COLORECT: try subsetting to specific mut type (along lines of SBS10a: REF = C vs SBS28: REF = T): 10a C>A, 10b C>T, 28 T>G
melahitsbak <- melahits
# melahits <- melahitsbak
# melahits <- melahits[melahits$ref == "C" & melahits$alt == "T"]


# freqtab <- c(table(factor(melahits$mut, levels = unique(samplesnvs$mut)))/table(factor(samplesnvs$mut, levels = unique(samplesnvs$mut))))
# samplefreqtab <- c(table(factor(melahits$sampleid, levels = unique(samplesnvs$sampleid)))/table(factor(samplesnvs$sampleid, levels = unique(samplesnvs$sampleid))))
# normlevels <- unique(paste0(samplesnvs$mut, "_", samplesnvs$sampleid))
# normfact <- c(table(factor(paste0(melahits$mut, "_", melahits$sampleid), levels = normlevels))/table(factor(paste0(samplesnvs$mut, "_", samplesnvs$sampleid), levels = normlevels)))
# # freqtab <- freqtab/sum(freqtab)
# sampledsnvs <- samplesnvs[sample(x = 1:length(samplesnvs), size = 10*length(melahits), replace = T, prob = normfact[paste0(samplesnvs$mut, "_", samplesnvs$sampleid)]), ]

ctrlseq <- samplesnvs$heptanuc #[samplesnvs$REF == "C"]
names(ctrlseq) <- paste0("ctrl", 1:length(ctrlseq))
writeXStringSet(x = ctrlseq, filepath = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/", ttype, "new_ctrl_seqs_over10.fasta"), format="fasta")
# writeXStringSet(x = ctrlseq, filepath = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/colorect_ctrl_seqs.fasta", format="fasta")
# writeXStringSet(x = ctrlseq, filepath = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/colorect_ctrl_seqs_SBS10b.fasta", format="fasta")

melaseq <- melahits$heptanuc #[melahits$ref == "C"]
names(melaseq) <- paste0("motif", 1:length(melaseq))
writeXStringSet(x = melaseq, filepath = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/", ttype, "_new_seqs_over10.fasta"), format="fasta")
# writeXStringSet(x = melaseq, filepath = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/colorect_seqs.fasta", format="fasta")
# writeXStringSet(x = melaseq, filepath = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/colorect_seqs_SBS10b.fasta", format="fasta")
}


### for CRC
# temp <- melahits$heptanuc[grep(pattern = "^[ACGT]{3}A[ACGT]TTC[GT]", x = melahits$heptanuc)]
# names(temp) <- paste0("motif", 1:length(temp))
# writeXStringSet(x = temp, filepath = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/colorect_seqs_motif1hits.fasta", format="fasta")

driverhits <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
driverhits[driverhits$sample %in% unique(melahits$sampleid) & grepl(pattern = "POLE", x = driverhits$gene_id), ]
driverhits[driverhits$sample %in% unique(melahits[grep(pattern = "CACACACA", x = melahits$heptanuc)]$sampleid), ]


### for CRC at TTCGA (run manual, need bits of code from above)
melahitsCpG <- melahits[grep(pattern = "^[ACGT]{5}TTCGA", x = melahits$heptanuc)]
sampledsnvsCpG <- sampledsnvs

# load bigwig with meth data
bw1 <- import.bw(con = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/28734.Roadmap.H-23769_729.WGB-Seq.signal.bigWig")
seqlevelsStyle(x = bw1) <- "Ensembl"
bw1 <- GRanges(bw1, seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
hitscores <- subsetByOverlaps(x = bw1, ranges = melahitsCpG)
hist(hitscores$score)
hist(bw1$score)
# sampledsnvsCpG <- 
sampledscores <- subsetByOverlaps(x = bw1, ranges = sampledsnvsCpG)
hist(sampledscores$score)

othermet <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/GSM1279519_CpGcontext.Colon.txt.gz", as.is = T, header = F)
bw1 <- GRanges(seqnames = othermet$V1, ranges = IRanges(start = othermet$V2, width = 1), mb = othermet$V3/(othermet$V3+othermet$V4), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))

hitscores <- subsetByOverlaps(x = bw1, ranges = melahitsCpG)
hist(hitscores$mb)
median(hitscores$mb, na.rm = T)
mean(hitscores$mb, na.rm = T)

sampledscores <- subsetByOverlaps(x = bw1, ranges = sampledsnvsCpG)
hist(sampledscores$mb)
median(sampledscores$mb, na.rm = T)
mean(sampledscores$mb, na.rm = T)

median(bw1$mb, na.rm = T)
mean(bw1$mb, na.rm = T)

wilcox.test(x = hitscores$mb, y = sampledscores$mb, alternative = "greater")
wilcox.test(x = hitscores$mb, y = bw1$mb[sample(x = 1:length(bw1), size = 100000, replace = F)], alternative = "greater")
hist(bw1$mb)





######## note, quick check for dihedrals + mid bond distances
params <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/structNFAT+DNA/distndihedrals.txt", as.is = T)
midbonds <- aggregate(x = params[, c("x", "y", "z")], by = list(molecule = params$molecule, pdbid = params$pdbid, chainid = params$chainid, res = params$resid), FUN = mean)
out <- by(data = midbonds[, c("x", "y", "z")], INDICES = paste0(midbonds$molecule, "_", midbonds$pdbid, "_", midbonds$chainid), FUN = function(x) sqrt(sum((x[1,] - x[2,])^2)), simplify = T)

outdf <- data.frame(do.call(rbind, strsplit(names(out), split = "_")), dist = c(out))
# outdf <- merge(x = outdf, y = params[, c("pdbid", "chainid", "dihedral")], by.x = c("X2", "X3"), by.y = c("pdbid", "chainid"), all.x = T, all.)
dihedrals <- params[!duplicated(paste0(params$pdbid, "_", params$chainid)), "dihedral"]
dihedrals <- setNames(object = dihedrals, nm = paste0(params$molecule, "_", params$pdbid, "_", params$chainid)[!duplicated(paste0(params$pdbid, "_", params$chainid))])
outdf$dihedral <- dihedrals[rownames(outdf)]


# df <- data.frame(dist = c(out), dih = dihedrals)[-length(out), ]

library(ggplot2)
p1 <- ggplot(data = outdf[outdf$X1 != "ideal", ], mapping = aes(y = dist, x = dihedral, colour = X1, group = X1, fill = X1)) 
p1 <- p1 + stat_ellipse(geom = "polygon", alpha = .1, show.legend = F) + labs(x = "dihedral angle", y = "distance")
p1 <- p1 + geom_jitter(width = .1, height = 0.025, alpha = .7, show.legend = F) + theme_minimal()
p1 <- p1 + scale_fill_manual(values = c(free = "#b2df8a", NFAT = "#a6cee3" ,ETS = "#1f78b4")) + scale_color_manual(values = c(free = "#b2df8a", NFAT = "#a6cee3" ,ETS = "#1f78b4"))
p1
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/dist_dih_plot_free-ETS-NFATc.pdf"), plot = p1, width = 6, height = 5, useDingbats=FALSE)
# sampledsnvs$trinuc <- getSeq(resize(allhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)


# parhitsmela <- parhits$heptanuc[parhits$histology_abbreviation == "ColoRect-AdenoCA"]
# names(parhitsmela) <- paste0("motif", 1:length(parhitsmela))
# writeXStringSet(x = parhitsmela, filepath = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/mela_seqs.fasta", format="fasta")
# 
# seqs <- paste0(">mot", 1:length(parhits), "\n", parhits$heptanuc)
# write(x = seqs, file = )

# 25 hits, 4 total from mnvs
paretshits <- parhits[grep(pattern = "[AT]TCC", x = parhits$heptanuc)]
any(!paretshits$sampleid %in% sampleids)
View(as.data.frame(paretshits))

table(paretshits$histology_abbreviation)
table(parhits$histology_abbreviation)

# parhits <- parhits[!seqnames(parhits) == "Y"]
# parhits <- parhits[parhits$sampleid %in% sampleids]

# 2 hits
thirdhits$trinuc <- getSeq(resize(thirdhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
thirdhits$heptanuc <- getSeq(resize(thirdhits, fix = "start", width = 7), x = BSgenome.Hsapiens.1000genomes.hs37d5)
thirdhits$histology_abbreviation <- sumtab$histology_abbreviation[match(x = thirdhits$sampleid, table = sumtab$samplename)]

allhits <- c(parhits[, c("ref", "alt", "trinuc")])


thirdetshits <- thirdhits[grep(pattern = "CTTCCG", x = thirdhits$heptanuc)]
any(!thirdetshits$sampleid %in% sampleids)
View(as.data.frame(thirdetshits))


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


driverhits <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
driverhits <- driverhits[driverhits$ref %in% c("A", "C", "G", "T") & driverhits$alt %in% c("A", "C", "G", "T"), c("chr", "pos")]
driverhits <- sort(unique(GRanges(seqnames = driverhits$chr, ranges = IRanges(start = driverhits$pos, width = 1))))
