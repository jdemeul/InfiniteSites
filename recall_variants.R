### check which PCAWG cohorts to recall variants from and dig out biallelics

# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(parallel)


RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"

read_pcawg_release_table <- function(release_table_file) {
  release_table <- read.delim(file = release_table_file, as.is = T)
  release_table <- release_table[release_table$wgs_exclusion_white_gray != "Excluded", ]
  splitaliquots <- strsplit(x = release_table$tumor_wgs_aliquot_id, split = ",")
  release_table_dedup <- release_table[rep(1:nrow(release_table), lengths(splitaliquots)), c("wgs_exclusion_white_gray", "normal_wgs_aliquot_id", "tumor_wgs_aliquot_id", "normal_wgs_bwa_alignment_bam_file_name", "tumor_wgs_bwa_alignment_bam_file_name", "dcc_project_code", "TiN")]
  release_table_dedup$tumor_wgs_aliquot_id <- unlist(splitaliquots)
  release_table_dedup$tumor_wgs_bwa_alignment_bam_file_name <- unlist(strsplit(x = release_table$tumor_wgs_bwa_alignment_bam_file_name, split = ","))
  return(release_table_dedup)
}


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

longdf <- expand.grid(trinucs$trinucleotides_mutations, trinucs$trinucleotides_mutations, stringsAsFactors = F)
colnames(longdf) <- c("mut1", "mut2")

longdf[, c("reftri1", "alt1")] <- do.call(rbind, strsplit(longdf$mut1, split = ">"))
longdf$ref1 <- substr(x = longdf$mut1, 2, 2)
longdf$alttri1 <- paste0(substr(x = longdf$mut1, 1, 1), longdf$alt1, substr(x = longdf$mut1, 3, 3))
longdf$alttri1 <- ifelse(longdf$alt1 %in% c("A", "G"), as.character(reverseComplement(DNAStringSet(longdf$alttri1))), longdf$alttri1)


longdf[, c("reftri2", "alt2")] <- do.call(rbind, strsplit(longdf$mut2, split = ">"))
longdf$ref2 <- substr(x = longdf$mut2, 2, 2)
longdf$alttri2 <- paste0(substr(x = longdf$mut2, 1, 1), longdf$alt2, substr(x = longdf$mut2, 3, 3))
longdf$alttri2 <- ifelse(longdf$alt1 %in% c("A", "G"), as.character(reverseComplement(DNAStringSet(longdf$alttri2))), longdf$alttri2)


longdf$mut1 <- factor(longdf$mut1, levels = trinucs$trinucleotides_mutations)
longdf$mut2 <- factor(longdf$mut2, levels = trinucs$trinucleotides_mutations)

#empty/impossible unless stated otherwise
longdf$type <- "empty"

# parallel defined as: 2x same reftri + alt
# third allele: 2x same reftri + different alt

# forward: alttri = reftri & alt2 != ref1
# back: alttri = reftri & alt2 = ref1

longdf$type <- ifelse(longdf$reftri1 == longdf$reftri2,
                      ifelse(longdf$alt1 == longdf$alt2, "parallel", "third_allele"),
                      ifelse(longdf$alttri1 == longdf$reftri2, ifelse(longdf$alttri2 == longdf$reftri1, "back", "forward"), longdf$type))
longdf$type <- factor(longdf$type, levels = c("back", "forward", "parallel", "third_allele", "empty"))
# longdf <- longdf[!(longdf$type == "third_allele" & as.integer(longdf$mut1) > as.integer(longdf$mut2)), ]

temp <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown/0009b464-b376-4fbc-8a56-da538269a02f_infsites_backfwd_allelic.txt", as.is = T)
all(as.character(longdf$mut1) == temp$mut1[-9217])
all(as.character(longdf$mut2) == temp$mut2[-9217])

third_alelle_idxs <- which(longdf$type == "third_allele")
forward_idxs <- which(longdf$type == "forward")





infsitesfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown/", pattern = "_infsites_backfwd_allelic.txt", full.names = T)
infsitesnums <- do.call(rbind, mclapply(X = infsitesfiles, FUN = function(x, taidx, fwdidx) {
  y <- read.delim(x, as.is = T)
  c(fwd = sum(y[fwdidx, "freqbf_hi"]), thirdallele = sum(y[taidx, "freqal_hi"]))
  }, taidx = third_alelle_idxs, fwdidx = forward_idxs, mc.preschedule = T, mc.cores = 16))
infsitesnums <- as.data.frame(infsitesnums)
infsitesnums$sampleid <- gsub(pattern = "_infsites_backfwd_allelic.txt", replacement = "", x = basename(infsitesfiles))
colnames(infsitesnums) <- c("mut1", "mut2", "freqbf_low", "freqbf_med", "freqbf_hi", "freqal_low", "freqal_med", "freqal_hi", "sampleid")

rslurmdf <- read_pcawg_release_table(release_table_file = RELEASETABLEFILE)
hitters <- infsitesnums[infsitesnums$fwd > 0 | infsitesnums$thirdallele > 0, "sampleid"]
outtab <- table(rslurmdf$dcc_project_code, rslurmdf$tumor_wgs_aliquot_id %in% hitters)
sort(outtab[,2]/rowSums(outtab))

# three top cohorts are COAD-US, SKCM-US, MELA-AU
rslurmdf <- rslurmdf[(rslurmdf$dcc_project_code %in% c("MELA-AU") | rslurmdf$tumor_wgs_aliquot_id %in% hitters) &
                       rslurmdf$wgs_exclusion_white_gray == "Whitelist" & (is.na(rslurmdf$TiN) | rslurmdf$TiN == 0) ,  ]

# prepare all the normals
# any(duplicated(rslurmdf$normal_wgs_bwa_alignment_bam_file_name))

write.table(x = rslurmdf[, 2:6], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/samples_to_recall.txt", quote = F, sep = "\t", row.names = F, col.names = T)

