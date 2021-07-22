

load_reference_alleles <- function(refallelesdir = REFALLELESDIR, chrominfo = genome_seqinfo) {
  require(readr)
  require(GenomicRanges)
  
  refallelefiles <- paste0(file.path(refallelesdir, "1000genomesAlleles2012_chr"), 1:23, ".txt")
  refalleles <- lapply(X = refallelefiles, function(x) readr::read_tsv(file = x, col_types = "iii"))
  chr <- rep(c(1:22, "X"), sapply(X = refalleles, FUN = nrow))
  refalleles <- as.data.frame(do.call(what = rbind, args = refalleles))
  refalleles_gr <- GenomicRanges::GRanges(seqnames = chr, IRanges::IRanges(start = refalleles$position, end = refalleles$position), mcols = refalleles[, -1], seqinfo = chrominfo)
  colnames(S4Vectors::mcols(refalleles_gr)) <- c("ref", "alt")
  return(refalleles_gr)
}


load_problematic_loci <- function(problocifile = PROBLOCIFILE, chrominfo = genome_seqinfo) {
  problemSNPs <- readr::read_tsv(file = problocifile, col_types = "ci")
  problemSNPs_gr <- GenomicRanges::GRanges(seqnames = problemSNPs$Chr, ranges = IRanges::IRanges(problemSNPs$Pos, end = problemSNPs$Pos), seqinfo = chrominfo)
  return(problemSNPs_gr)
}

# 
# 
# load_allelecounts <- function(allelecountsfile, chrominfo = genome_seqinfo) {
#   require(readr)
#   require(GenomicRanges)
#   
#   allelecounts <- readr::read_tsv(file = allelecountsfile, col_types = "ciiiii")
#   allelecounts_gr <- GenomicRanges::GRanges(seqnames = allelecounts$Chromosome, IRanges::IRanges(start = allelecounts$Position, end = allelecounts$Position), mcols = allelecounts[, -c(1,2)], seqinfo = chrominfo)
#   colnames(S4Vectors::mcols(allelecounts_gr)) <- c("mutCountT1", "mutCountT2", "mutCountN1", "mutCountN2")
#   
#   return(allelecounts_gr)
# }


load_allelecounts <- function(allelecountsfile, scratchdir) {
  require(readr)
  
  # avoid interference by same samples
  scratchdir_bc <- file.path(scratchdir, paste0(sample(x = c(letters), size = 12, replace = T), collapse = ""))
  dir.create(path = scratchdir_bc, recursive = T)
  
  # unpack allelecounts
  system(command = paste0("tar -xzf ", allelecountsfile, " -C ", scratchdir_bc), wait = T)
  
  # read allelecounts
  sampleid <- sub(pattern = "_allelecounts.tar.gz", replacement = "", x = basename(allelecountsfile))
  allelecountsfolder <- file.path(scratchdir_bc, paste0(sampleid, "_allelecounts"))
  
  allelecounts <- do.call(rbind, lapply(X = file.path(allelecountsfolder, paste0(sampleid, "_alleleFrequencies_chr", 1:23, ".txt")),
                                        function(x) readr::read_tsv(file = x, col_types = "ciiiiii", col_names = c("chr", "pos", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth"), comment = "#")))
  
  # clean up
  unlink(x = scratchdir_bc, recursive = T)
  
  return(allelecounts)
}


load_allelecounts_normal <- function(allelecountsfile_normal, alleles, chrominfo = genome_seqinfo, actempdir = ALLELECOUNTSTEMPDIR) {
  require(readr)
  require(GenomicRanges)

  allelecounts <- load_allelecounts(allelecountsfile = allelecountsfile_normal, scratchdir = actempdir)

  allelecounts_gr <- GenomicRanges::GRanges(seqnames = allelecounts$chr, IRanges::IRanges(start = allelecounts$pos, end = allelecounts$pos), seqinfo = chrominfo)
  mcols(allelecounts_gr) = allelecounts[, -c(1,2)]
  allelecounts_gr <- sort(allelecounts_gr)

  # allelehits <- findOverlaps(query = allelecounts_gr, subject = alleles)
  # allelecounts_gr <- allelecounts_gr[subjectHits(allelehits)]
  mcols(allelecounts_gr)$ref <- factor(x = mcols(alleles)$ref, levels = 1:4, labels = c("A", "C", "G", "T"))
  mcols(allelecounts_gr)$alt <- factor(x = mcols(alleles)$alt, levels = 1:4, labels = c("A", "C", "G", "T"))
  mcols(allelecounts_gr)$mutCountN1 <- ifelse(mcols(allelecounts_gr)$ref == "A", mcols(allelecounts_gr)$Count_A,
                                              ifelse(mcols(allelecounts_gr)$ref == "C", mcols(allelecounts_gr)$Count_C,
                                                     ifelse(mcols(allelecounts_gr)$ref == "G", mcols(allelecounts_gr)$Count_G,
                                                            mcols(allelecounts_gr)$Count_T)))
  mcols(allelecounts_gr)$mutCountN2 <- ifelse(mcols(allelecounts_gr)$alt == "A", mcols(allelecounts_gr)$Count_A,
                                              ifelse(mcols(allelecounts_gr)$alt == "C", mcols(allelecounts_gr)$Count_C,
                                                     ifelse(mcols(allelecounts_gr)$alt == "G", mcols(allelecounts_gr)$Count_G,
                                                            mcols(allelecounts_gr)$Count_T)))
  mcols(allelecounts_gr) <- mcols(allelecounts_gr)[, c("mutCountN1", "mutCountN2", "ref", "alt")]
  
  # colnames(S4Vectors::mcols(allelecounts_gr)) <- c("mutCountT1", "mutCountT2", "mutCountN1", "mutCountN2")
  
  return(allelecounts_gr)
}



get_clean_hetSNPs <- function(sampleid, allelecountsfile, indelfile, breaksfile, probloci = NULL,
                                   chrominfo = genome_seqinfo, HDI90 = c(lower = 1/3, upper = 2/3), mindepth = 10, alleles,
                                   outfile, indelfilter = T, breakfilter = T, hdifilter = T, heterozygousfilter = .1) {
  require(GenomicRanges)
  require(readr)

  # allelecounts_gr <- load_allelecounts(allelecountsfile = allelecountsfile, chrominfo = chrominfo)
  allelecounts_gr <- load_allelecounts_normal(allelecountsfile = allelecountsfile, alleles = alleles, chrominfo = chrominfo)
  
  # basic subset by coverage of both alleles, to quickly shrink data down
  normalcov <- mcols(allelecounts_gr)$mutCountN1 + mcols(allelecounts_gr)$mutCountN2
  normalbaf <- mcols(allelecounts_gr)$mutCountN1/normalcov
  allelecounts_gr <- allelecounts_gr[which(normalcov >= mindepth & normalbaf >= heterozygousfilter & normalbaf <= 1-heterozygousfilter) ]
  # allelecounts_gr <- allelecounts_gr[mcols(allelecounts_gr)$mutCountN1 >= minalleledepth & mcols(allelecounts_gr)$mutCountN2 >= minalleledepth]
  
  # Remove general problematic loci
  if (!is.null(probloci)) {
    isbadSNP <- IRanges::overlapsAny(query = allelecounts_gr, subject = probloci)
  } else {
    isbadSNP <- F
  }
  
  # Remove SNPs too close to indels
  # indelfile <- list.files(path = indeldir, pattern = paste0(sampleid, ".consensus.20161006.somatic.indel.vcf.gz$"), full.names = T, recursive = T)
  if (length(indelfile) == 1 & indelfilter) {
    indels <- SummarizedExperiment::rowRanges(VariantAnnotation::readVcf(file = indelfile), genome = chrominfo)
    isnearindel <- IRanges::overlapsAny(query = allelecounts_gr, subject = indels, maxgap = 100)
  } else {
    isnearindel <- F
  }
  
  # Remove SNPs too close to SV/CN breaks
  # breaksfile <- list.files(path = breakpointsdir, pattern = paste0(sampleid, ".consensus.20170119.somatic.cna.annotated.txt$"), full.names = T)
  if (length(breaksfile) == 1 & breakfilter) {
    breaks <- read.delim(file = breaksfile, as.is = T)
    breaks_gr <- c(GRanges(seqnames = breaks$chromosome, ranges = IRanges(start = breaks$start, end = breaks$start), seqinfo = chrominfo), 
                  GRanges(seqnames = breaks$chromosome, ranges = IRanges(start = breaks$end, end = breaks$end), seqinfo = chrominfo))
    isnearbreak <- IRanges::overlapsAny(query = allelecounts_gr, subject = breaks_gr, maxgap = 100)
  } else {
    isnearbreak <- F
  }
  
  allelecounts_gr <- allelecounts_gr[!(isbadSNP | isnearindel | isnearbreak)]
  

  # only retain clean, confident het SNPs
  if (hdifilter) {
    pconfint <- mapply(FUN = function(n1, n2) qbeta(c(0.05, 0.95), shape1 = n1+2, shape2 = n2+2, lower.tail = TRUE, log.p = FALSE), n1 = mcols(allelecounts_gr)$mutCountN1, n2 = mcols(allelecounts_gr)$mutCountN2)
    allelecounts_gr <- allelecounts_gr[pconfint[1,] >= HDI90[["lower"]] & pconfint[1,] < 0.5 & pconfint[2,] <= HDI90[["upper"]] & pconfint[2,] > 0.5]
  }
  
  # ## add in ref and alt alleles for the remaining SNPs
  # annotHits <- findOverlaps(query = allelecounts_gr, subject = alleles)
  # # check whether all used alleles are indeed present:
  # # all(queryHits(annotHits) %in% 1:length(allelecounts_gr))
  # mcols(allelecounts_gr)[, c("ref", "alt")] <- mcols(alleles)[subjectHits(annotHits), ]
  # mcols(allelecounts_gr)$ref <- factor(mcols(allelecounts_gr)$ref, levels = 1:4, labels = c("A", "C", "G", "T"))
  # mcols(allelecounts_gr)$alt <- factor(mcols(allelecounts_gr)$alt, levels = 1:4, labels = c("A", "C", "G", "T"))
  
  ## writing out loci file of good SNPs for use in phasing code
  hetSNPs <- as.data.frame(allelecounts_gr)[ , c("seqnames", "start", "ref", "alt")]
  colnames(hetSNPs) <- c("chr", "pos", "ref", "alt")
  readr::write_tsv(x = hetSNPs, path = outfile, col_names = T)
  
  return(NULL)
}




