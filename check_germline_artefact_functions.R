#### additional checks to exclude contamination by germline calls
# 


#### function to annotate PCAWG vcf with the additional QC metrics (carry this along to the end)
## get allele counts in matched normal

annot_normal_read_counts <- function(normalbam, snvs_vcf, minbq = 20, minmq = 35, ncores = NCORES) {
  
  if (!file.exists(normalbam)) {
    print("Normal bam not found")
    snvs_vcf$n_ref_count <- NA
    snvs_vcf$n_alt_count <- NA
    snvs_vcf$n_total_cov <- NA
    return(snvs_vcf)
  }
  
  # pile up parameters
  flag <- scanBamFlag(isNotPassingQualityControls = F, isDuplicate = F, isProperPair = T)
  # scbamparam <- ScanBamParam(which = snvs_vcf, flag = flag) #, mapqFilter=minMapQ # simpleCIGAR: allow only reads with M in CIGAR string, no clipping, indels or whatever
  
  pupparam <- PileupParam(max_depth=250, min_base_quality=minbq, min_mapq=minmq, min_nucleotide_depth=0,
                          min_minor_allele_depth=0, distinguish_strands=F, distinguish_nucleotides=TRUE,
                          ignore_query_Ns=TRUE, include_deletions=F, include_insertions=F, left_bins=NULL, query_bins=NULL, cycle_bins=NULL)
  # actual pile up -- leaner versions
  vcflocils <- split(x = granges(snvs_vcf), f = seqnames(snvs_vcf))
  vcflocils <- vcflocils[sapply(X = vcflocils, FUN = length) > 0]
  
  if (ncores == 1 || length(snvs_vcf) < 1000) {
    pupout <- lapply(X = vcflocils, FUN = function(vcfpos, flag, normalbam, pupparam) {
      # browser()
      scbamparam <- ScanBamParam(which = vcfpos, flag = flag) #, mapqFilter=minMapQ # simpleCIGAR: allow only reads with M in CIGAR string, no clipping, indels or whatever
      pupout <- pileup(file = normalbam, index = normalbam, scanBamParam = scbamparam, pileupParam = pupparam)
      return(pupout)
    }, flag = flag, normalbam = normalbam, pupparam = pupparam)
  } else {
      pupout <- mclapply(X = vcflocils, FUN = function(vcfpos, flag, normalbam, pupparam) {
      # browser()
      scbamparam <- ScanBamParam(which = vcfpos, flag = flag) #, mapqFilter=minMapQ # simpleCIGAR: allow only reads with M in CIGAR string, no clipping, indels or whatever
      pupout <- pileup(file = normalbam, index = normalbam, scanBamParam = scbamparam, pileupParam = pupparam)
      return(pupout)
    }, flag = flag, normalbam = normalbam, pupparam = pupparam, mc.cores = ceiling(ncores)/2, mc.preschedule = F)
  }
  

  # pupout <- pileup(file = normalbam, index = normalbam, scanBamParam = scbamparam, pileupParam = pupparam)
  # reshape output
  pupout <- dcast(data = do.call(rbind, pupout), formula = seqnames + pos ~ nucleotide, value.var = "count", fill = 0)
  pupout_gr <- GRanges(seqnames = pupout$seqnames, ranges = IRanges(pupout$pos, width = 1))
  sharedcols <- intersect(x = c("A", "C", "G", "T"), colnames(pupout))
  mcols(pupout_gr)[, c("A", "C", "G", "T")] <- 0
  mcols(pupout_gr)[, sharedcols] <- pupout[, sharedcols]
  if (nrow(pupout) == 1) {
    pupout_gr$total <- sum(pupout[, -c(1,2)])
  } else if (length(sharedcols) == 1) {
    pupout_gr$total <- pupout[, sharedcols]
  } else {
    pupout_gr$total <- rowSums(pupout[, -c(1,2)])
  }

  
  # enrich snvs_gr
  matchidxs <- nearest(x = snvs_vcf, subject = pupout_gr)
  
  snvs_vcf$n_ref_count <- 0
  snvs_vcf$n_alt_count <- 0
  snvs_vcf$n_total_cov <- 0
  
  # had an NA here, no read counts collected at locus
  goodidxs <- which(!is.na(matchidxs))
  matchidxs <- na.omit(matchidxs)
  
  # modified to take in biallelic divergent variants
  snvs_vcf$n_ref_count[goodidxs] <- ifelse(snvs_vcf[goodidxs]$REF == "A", pupout_gr[matchidxs]$A, 
                                           ifelse(snvs_vcf[goodidxs]$REF == "C", pupout_gr[matchidxs]$C,
                                                  ifelse(snvs_vcf[goodidxs]$REF == "G", pupout_gr[matchidxs]$G, pupout_gr[matchidxs]$'T')))
  
  if (any(lengths(snvs_vcf[goodidxs]$ALT) > 1)) {
    alt1 <- unlist(DNAStringSetList(sapply(snvs_vcf[goodidxs]$ALT, "[", 1)))
    alt2 <- unlist(DNAStringSetList(sapply(snvs_vcf[goodidxs]$ALT, "[", 2)))
    snvs_vcf$n_alt_count[goodidxs] <- ifelse(alt1 == "A", pupout_gr[matchidxs]$A, 
                                             ifelse(alt1 == "C", pupout_gr[matchidxs]$C,
                                                    ifelse(alt1 == "G", pupout_gr[matchidxs]$G, pupout_gr[matchidxs]$'T')))
    
    snvs_vcf$n_alt_count2 <- 0
    snvs_vcf$n_alt_count2[goodidxs] <- ifelse(alt2 == "A", pupout_gr[matchidxs]$A, 
                                             ifelse(alt2 == "C", pupout_gr[matchidxs]$C,
                                                    ifelse(alt2 == "G", pupout_gr[matchidxs]$G, pupout_gr[matchidxs]$'T')))
  } else {
    snvs_vcf$n_alt_count[goodidxs] <- ifelse(snvs_vcf[goodidxs]$ALT == "A", pupout_gr[matchidxs]$A, 
                                             ifelse(snvs_vcf[goodidxs]$ALT == "C", pupout_gr[matchidxs]$C,
                                                    ifelse(snvs_vcf[goodidxs]$ALT == "G", pupout_gr[matchidxs]$G, pupout_gr[matchidxs]$'T')))
  }
  
  snvs_vcf$n_total_cov[goodidxs] <- pupout_gr[matchidxs]$total
  
  return(snvs_vcf)
}


# already annotated with 1000G AF even
# annot_1kG_matches <- function(refalleles_gr, snvs_vcf) {
#   
#   refoverlaps <- findOverlaps(query = snvs_vcf, subject = refalleles_gr, type = "equal")
#   snvs_vcf$in1kG <- F
#   snvs_vcf[queryHits(refoverlaps)]$in1kG <- snvs_vcf[queryHits(refoverlaps)]$ALT == refalleles_gr[subjectHits(refoverlaps)]$alt
#   
#   return(snvs_vcf)
# }


annot_hg38_consistency <- function(snvs_vcf, hg19tohg38) {
  
  allsnvs_gr <- snvs_vcf
  seqlevelsStyle(allsnvs_gr) <- "UCSC"
  allsnvs_grhg38 <- liftOver(chain = hg19tohg38chain, x = allsnvs_gr)
  # seqlevelsStyle(allsnvs) <- "ensembl"
  allsnvs_grhg38 <- unlist(allsnvs_grhg38[names(which(lengths(allsnvs_grhg38) == 1))], use.names = T)
  names(allsnvs_grhg38) <- unlist(sapply(X = strsplit(x = names(allsnvs_grhg38), split = ".", fixed = T), FUN = "[", 1)) # fix names if appended - package version differences
  seqlevelsStyle(allsnvs_grhg38) <- "ensembl"
  allsnvs_grhg38$REF38 <- getSeq(x = BSgenome.Hsapiens.NCBI.GRCh38, names = allsnvs_grhg38)
  # goodloci <- names(allsnvshg38)[which(allsnvshg38$REF38 == allsnvshg38$REF)]
  
  # check consistency in GRCh38
  snvs_vcf$hg38clean <- F
  refmatchidxs <- which(allsnvs_grhg38$REF38 == allsnvs_grhg38$REF)
  if (length(refmatchidxs) > 0) {
    # snvs_vcf[names(allsnvs_grhg38)[refmatchidxs]]$hg38clean <- T
    mcols(snvs_vcf)[names(allsnvs_grhg38)[refmatchidxs], "hg38clean"] <- T
  }
  
  return(snvs_vcf)
  
}

# if "Validation_status field available, USE!!!
annot_validation_status <- function(snvs_vcf) {
  
  if (!"Validation_status" %in% colnames(mcols(snvs_vcf))) {
    snvs_vcf$Validation_status <- NA
  }

  return(snvs_vcf)
}



library(rtracklayer)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(reshape2)


add_snv_qc <- function(sampleid, sampledir, releasetable, snvs_vcf, hg19tohg38, normalbam, minbq = 20, minmq = 35, ncores = NCORES, bamdir = BAMDIR, checkbam = T) {
  
  annotfile <- file.path(sampledir, paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.annot.gr.RDS"))
  
  if (file.exists(annotfile)) {
    snvs_vcf <- readRDS(file = annotfile)
  } else {
    if (checkbam) {
      # annotate matched normal read counts/alt alleles
      if (grepl(pattern = "/srv/shared/vanloo/ICGC/", x = bamdir, fixed = T)) {
        normalbam <- file.path(bamdir, paste0(releasetable[releasetable$tumor_wgs_aliquot_id == sampleid, c("dcc_project_code", "normal_wgs_bwa_alignment_bam_file_name")], collapse = "/WGS/"))
      } else {
        normalbam <- file.path(bamdir, releasetable[releasetable$tumor_wgs_aliquot_id == sampleid, "dcc_project_code"], sampleid, paste0(sampleid, "_normal.aln.recal.bam"))
      }
      snvs_vcf <- annot_normal_read_counts(normalbam = normalbam, snvs_vcf = snvs_vcf, minbq = minbq, minmq = minmq, ncores = ncores)
    } else {
      snvs_vcf$n_ref_count <- NA
      snvs_vcf$n_alt_count <- NA
      snvs_vcf$n_total_cov <- NA
    }

    # anntoate validation status
    snvs_vcf <- annot_validation_status(snvs_vcf = snvs_vcf)
    
    # annotate hg38 consistency
    snvs_vcf <- annot_hg38_consistency(snvs_vcf = snvs_vcf, hg19tohg38 = hg19tohg38)
    
    # check for presence of QC columns
    if (!all(c("n_ref_count", "n_alt_count", "n_total_cov",
               "hg38clean", "Validation_status") %in% colnames(mcols(snvs_vcf)))) {
      print("Some QC columns are missing")
    }
    saveRDS(object = snvs_vcf, file = annotfile)
  }
  
  return(snvs_vcf)
}






