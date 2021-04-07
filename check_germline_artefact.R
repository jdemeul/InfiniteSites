#### additional checks to exclude contamination by germline calls
# 
# library(VariantAnnotation)
# library(BSgenome.Hsapiens.1000genomes.hs37d5)
# library(BSgenome.Hsapiens.NCBI.GRCh38)
# library(ggplot2)
# library(Rsamtools)
# library(reshape2)
# 
# 
# SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/"
# SNVMNVINDEANNOTLDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/normaldepthannotatedvcfs/"
# RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
# SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
# PARHITSFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt"
# DIVHITSFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt"
# hg19tohg38chainfile <- "/srv/shared/vanloo/home/jdemeul/projects/2019_Rambow_MELA/data/hg19ToHg38.over.chain"
# 
# 
# sumtab <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
# releasetab <- read.delim(file = RELEASETABLEFILE, as.is = T)
# samplesofinterest <- unique(c(read.delim(file = PARHITSFILE, as.is = T)$sampleid, read.delim(file = DIVHITSFILE, as.is = T)$sampleid))
# parhits <- read.delim(file = PARHITSFILE, as.is = T)
# divhits <- read.delim(file = DIVHITSFILE, as.is = T)
# hg19tohg38chain <- import.chain(con = hg19tohg38chainfile)
# 
# parhits$id <- paste0(parhits$chr, ":", parhits$start, "_", parhits$ref, "/", parhits$alt)
# divhits$id <- paste0(divhits$seqnames, ":", divhits$start, "_", divhits$REF, "/", divhits$in_PCAWG)
# 
# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# # sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
# 
# # check depth of coverage in the normal (get info from vcf)
# get_snvs <- function(sampleid, snvsmnvdir) {
#   # snv_mnvfile <- file.path(snvsmnvdir, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz"))
#   snv_mnvfile <- file.path(snvsmnvdir, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.wdepth.vcf"))
#   snvs <- readVcf(file = snv_mnvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
#   
#   return(snvs)
# }
# 
# allsnvs <- get_snvs(sampleid = sampleid, snvsmnvdir = SNVMNVINDEANNOTLDIR)
# 
# allsnvs_gr <- granges(allsnvs)
# allsnvs_gr$ndepth <- info(allsnvs)$BAM_DEPTH
# allsnvs_gr$nearindel <- info(allsnvs)$snv_near_indel
# seqlevelsStyle(allsnvs_gr) <- "UCSC"
# allsnvs_grhg38 <- liftOver(chain = hg19tohg38chain, x = allsnvs_gr)
# # seqlevelsStyle(allsnvs) <- "ensembl"
# allsnvs_grhg38 <- unlist(allsnvs_grhg38[which(lengths(allsnvs_grhg38) == 1)], use.names = F)
# seqlevelsStyle(allsnvs_grhg38) <- "ensembl"
# allsnvs_grhg38$REF38 <- getSeq(x = BSgenome.Hsapiens.NCBI.GRCh38, names = allsnvs_grhg38)
# # goodloci <- names(allsnvshg38)[which(allsnvshg38$REF38 == allsnvshg38$REF)]
# 
# # check consistency in GRCh38
# allsnvs_gr$hg38clean <- F
# allsnvs_gr[names(allsnvs_grhg38)[which(allsnvs_grhg38$REF38 == allsnvs_grhg38$REF)]]$hg38clean <- T
# allsnvs_gr$parallel <- F
# allsnvs_gr[parhits[parhits$sampleid == sampleid, "id"]]$parallel <- T
# 
# table(biallelic = allsnvs_gr$parallel, hg38clean = allsnvs_gr$hg38clean)
# chisq.test(table(biallelic = allsnvs_gr$parallel, hg38clean = allsnvs_gr$hg38clean))
# ### conclusion, remove the few h38badloci
# ### bump this bit of code into call_parallel_violations (GC_by_AF functions)
# ###
# 
# allsnvs_gr$divergent <- F
# allsnvs_gr[divhits[divhits$sampleid == sampleid & divhits$in_PCAWG != "", "id"]]$divergent <- T
# 
# # check enrichment for known SNPs
# allsnvs_gr$dbsnp <- F
# allsnvs_gr[names(allsnvs)[!is.na(info(allsnvs)$dbsnp) & !info(allsnvs)$dbsnp_somatic]]$dbsnp <- T
# 
# table(biallelic = allsnvs_gr$parallel, dbsnp = allsnvs_gr$dbsnp)
# chisq.test(table(biallelic = allsnvs_gr$parallel, dbsnp = allsnvs_gr$dbsnp), simulate.p.value = T, B = 1e6)
# 
# # check whether near indel
# table(biallelic = allsnvs_gr$parallel, nearindel = allsnvs_gr$nearindel)
# chisq.test(table(biallelic = allsnvs_gr$parallel, nearindel = allsnvs_gr$nearindel))
# 
# # check normal depth
# allsnvs_df <- as.data.frame(mcols(allsnvs_gr)[, c("ndepth", "parallel", "divergent")])
# p1 <- ggplot(data = allsnvs_df, mapping = aes(x = ndepth)) + geom_histogram(mapping = aes(fill = parallel), alpha = .6, position = "identity") + scale_y_log10() + lims(x = c(0,60))
# p1

# pbinom(q = 0, size = 15, prob = .5, lower.tail = T)

# check in PCAWG pilot validation results for deep-seq

# check in largest cohorts (MELA and COAD whether any differences)



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
hg19tohg38chainfile <- "/srv/shared/vanloo/home/jdemeul/projects/2019_Rambow_MELA/data/hg19ToHg38.over.chain"
hg19tohg38chain <- import.chain(con = hg19tohg38chainfile)


add_snv_qc <- function(sampleid, sampledir, releasetable, snvs_vcf, hg19tohg38, normalbam, minbq = 20, minmq = 35, ncores = NCORES, bamdir = BAMDIR) {
  
  annotfile <- file.path(sampledir, paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.annot.gr.RDS"))
  
  if (file.exists(annotfile)) {
    snvs_vcf <- readRDS(file = annotfile)
  } else {
    # annotate matched normal read counts/alt alleles
    if (grepl(pattern = "/srv/shared/vanloo/ICGC/", x = bamdir, fixed = T)) {
      normalbam <- file.path(bamdir, paste0(releasetable[releasetable$tumor_wgs_aliquot_id == sampleid, c("dcc_project_code", "normal_wgs_bwa_alignment_bam_file_name")], collapse = "/WGS/"))
    } else {
      normalbam <- file.path(bamdir, releasetable[releasetable$tumor_wgs_aliquot_id == sampleid, "dcc_project_code"], sampleid, paste0(sampleid, "_normal.aln.recal.bam"))
    }
    snvs_vcf <- annot_normal_read_counts(normalbam = normalbam, snvs_vcf = snvs_vcf, minbq = minbq, minmq = minmq, ncores = ncores)
    
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






