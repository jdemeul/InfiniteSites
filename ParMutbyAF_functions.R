### functions for identifying gene conversion based on the allele frequencies

# get tumour allelecounts 
get_tum_allecounts_gr <- function(allelecountsdir, sampleid, bsgenome, reference_alleles, tempdir = TEMPDIR) {
  allelecountsfile_tum <- file.path(allelecountsdir, paste0(sampleid, "_beagle5"), paste0(sampleid, "_alleleFrequencies_chr", c(1:23), ".txt"))
  if (all(file.exists(allelecountsfile_tum))) {
    # allelecounts_tum <- load_allelecounts_notarball(allelecountsfile = allelecountsfile_tum, scratchdir = tempdir)
    allelecounts_tum <- do.call(rbind, lapply(X = allelecountsfile_tum,
                                          function(x) readr::read_tsv(file = x, col_types = "ciiiiii", col_names = c("chr", "pos", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth"), comment = "#")))
  } else {
    return(NULL)
  }
  allelecounts_tum_gr <- GRanges(seqnames = allelecounts_tum$chr, IRanges(start = allelecounts_tum$pos, end = allelecounts_tum$pos), seqinfo = seqinfo(bsgenome))
  mcols(allelecounts_tum_gr) <- allelecounts_tum[, -c(1,2)]
  
  allelematches <- findOverlaps(query = allelecounts_tum_gr, subject = reference_alleles)
  allelecounts_tum_gr <- allelecounts_tum_gr[queryHits(allelematches)]
  mcols(allelecounts_tum_gr)[, c("ref", "alt")] <- mcols(reference_alleles[subjectHits(allelematches)])

  
  mcols(allelecounts_tum_gr)$refCount <- ifelse(mcols(allelecounts_tum_gr)$ref == "A", mcols(allelecounts_tum_gr)$Count_A,
                                                ifelse(mcols(allelecounts_tum_gr)$ref == "C", mcols(allelecounts_tum_gr)$Count_C,
                                                       ifelse(mcols(allelecounts_tum_gr)$ref == "G", mcols(allelecounts_tum_gr)$Count_G,
                                                              mcols(allelecounts_tum_gr)$Count_T)))
  mcols(allelecounts_tum_gr)$altCount <- ifelse(mcols(allelecounts_tum_gr)$alt == "A", mcols(allelecounts_tum_gr)$Count_A,
                                                ifelse(mcols(allelecounts_tum_gr)$alt == "C", mcols(allelecounts_tum_gr)$Count_C,
                                                       ifelse(mcols(allelecounts_tum_gr)$alt == "G", mcols(allelecounts_tum_gr)$Count_G,
                                                              mcols(allelecounts_tum_gr)$Count_T)))
  
  mcols(allelecounts_tum_gr) <- mcols(allelecounts_tum_gr)[, -c(1:5)]
  return(allelecounts_tum_gr)
}



load_1000G_reference_alleles_new <- function(refallelesdir, chrominfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)) {
  # refallelefiles <- paste0(refallelesdir, "1kg.phase3.v5a_GRCh37nounref_allele_index_chr", c(1:22, "X"), ".txt")
  refallelefiles <- paste0(refallelesdir, "chr", c(1:22, "X"), ".1kg.phase3.v5a_GRCh37nounref_allele_index.txt")
  refalleles <- lapply(X = refallelefiles, function(x) read_tsv(file = x, col_types = "iii"))
  chr <- rep(c(1:22, "X"), sapply(X = refalleles, FUN = nrow))
  refalleles <- as.data.frame(do.call(what = rbind, args = refalleles))
  refalleles_gr <- GRanges(seqnames = chr, IRanges(start = refalleles$position, end = refalleles$position), seqinfo = chrominfo)
  mcols(refalleles_gr)$ref <- factor(refalleles$a0, levels = 1:4, labels = c("A", "C", "G", "T"))
  mcols(refalleles_gr)$alt <- factor(refalleles$a1, levels = 1:4, labels = c("A", "C", "G", "T"))
  return(refalleles_gr)
}


# get phased BAFs (for QC purposes, plotting, validation)
get_phased_BAF <- function(bafdir, sampleid, bsgenome, allelecounts) {
  phasedbaf_file <- file.path(bafdir, paste0(sampleid, "_beagle5"), paste0(sampleid, ".BAFsegmented.txt"))
  if (file.exists(phasedbaf_file)) {
    phasedbaf <- read_tsv(file = phasedbaf_file, col_types = "cinnn")
  } else {
    phasedbaf_file_alt <- file.path("/camp/project/proj-emedlab-vanloo/bafsegmented_JD", paste0(sampleid, ".BAFsegmented.txt.gz"))
    print(paste0("No segmented BAF data found at ", phasedbaf_file, ". \n Using segmented BAF data from ", phasedbaf_file_alt))
    if (file.exists(phasedbaf_file_alt)) {
      phasedbaf <- read_tsv(file = phasedbaf_file_alt, col_types = "cinnn")
    } else {
      return(NULL)
    }
  }
  phasedbaf_gr <- GRanges(seqnames = phasedbaf$Chromosome, IRanges(start = phasedbaf$Position, end = phasedbaf$Position), seqinfo = seqinfo(bsgenome))
  mcols(phasedbaf_gr) <- phasedbaf[, -c(1:2)]
  
  locimatches <- findOverlaps(query = phasedbaf_gr, subject = allelecounts)
  phasedbaf_gr <- phasedbaf_gr[queryHits(locimatches)]
  mcols(phasedbaf_gr)[, c("ref", "alt", "refCount", "altCount")] <- mcols(allelecounts)[subjectHits(locimatches), c("ref", "alt", "refCount", "altCount")]
  
  # use simple BAFphased to determine which is "major" allele
  mcols(phasedbaf_gr)$MajCount <- ifelse(mcols(phasedbaf_gr)$BAFphased > .5,
                                         pmax(mcols(phasedbaf_gr)$refCount, mcols(phasedbaf_gr)$altCount),
                                         pmin(mcols(phasedbaf_gr)$refCount, mcols(phasedbaf_gr)$altCount))
  mcols(phasedbaf_gr)$MinCount <- ifelse(mcols(phasedbaf_gr)$BAFphased > .5,
                                         pmin(mcols(phasedbaf_gr)$refCount, mcols(phasedbaf_gr)$altCount),
                                         pmax(mcols(phasedbaf_gr)$refCount, mcols(phasedbaf_gr)$altCount))
# 
#   mcols(phasedbaf_gr)$seg <- Rle(paste0(seqnames(phasedbaf_gr), "_", mcols(phasedbaf_gr)$BAFseg))
#   runValue(mcols(phasedbaf_gr)$seg) <- 1:nrun(mcols(phasedbaf_gr)$seg)
  
  return(phasedbaf_gr)
}


# get the GC corrected logR 
get_GCcorr_logr <- function(logrdir, sampleid, bsgenome) {
  gccorrlogr_file <- file.path(logrdir, paste0(sampleid, "_beagle5"), paste0(sampleid, "_mutantLogR_gcCorrected.tab"))
  if (file.exists(gccorrlogr_file)) {
    gccorrlogr <- read_tsv(file = gccorrlogr_file, col_types = "cin")
  } else {
    gccorrlogr_file_alt <- gsub(pattern = "_gcCorrected", replacement = "", fixed = T, x = gccorrlogr_file)
    if (file.exists(gccorrlogr_file_alt)) {
      gccorrlogr <- read_tsv(file = gccorrlogr_file_alt, col_types = "cin")
    } else {
      return(NULL)
    }
  }
  gccorrlogr_gr <- GRanges(seqnames = gccorrlogr$Chromosome, IRanges(start = gccorrlogr$Position, end = gccorrlogr$Position), logr = gccorrlogr[, 3, drop = T],  seqinfo = seqinfo(bsgenome))

  # # get logR points within the segments
  # logrhits <- findOverlaps(query = segments_gr, subject = gccorrlogr_gr)
  # gccorrlogr_gr <- gccorrlogr_gr[subjectHits(logrhits)] # not all are within BAF segments, male X has large part removed
  # mcols(gccorrlogr_gr)$seg <- queryHits(logrhits)
  
  return(gccorrlogr_gr)  
}


# summarise BAF/LogR per segment
summarise_baflogr_segments <- function(sampleid, sampledir, phasedbaf, logr, segments_gr, rhopsi) {

  # create baf/logr summary per segment
  mcols(segments_gr)$nhetsnps <- 0L
  mcols(segments_gr)$mu_logr <- as.numeric(rep(x = NA, length(segments_gr)))
  mcols(segments_gr)$rho_logr <- as.numeric(rep(x = NA, length(segments_gr)))
  mcols(segments_gr)$mu_baf <- as.numeric(rep(x = NA, length(segments_gr)))
  mcols(segments_gr)$rho_baf <- as.numeric(rep(x = NA, length(segments_gr)))
  
  logrsegoverlaps <- findOverlaps(query = segments_gr, subject = logr)
  
  mcols(segments_gr)[unique(queryHits(logrsegoverlaps)), "mu_logr"] <- by(data = mcols(logr)[subjectHits(logrsegoverlaps), "logr"], INDICES = queryHits(logrsegoverlaps), FUN = median, simplify = T, na.rm = T)
  mcols(segments_gr)[unique(queryHits(logrsegoverlaps)), "rho_logr"] <- by(data = mcols(logr)[subjectHits(logrsegoverlaps), "logr"], INDICES = queryHits(logrsegoverlaps), FUN = mad, simplify = T, na.rm = T)
  
  # mcols(segments_gr)$rho_logr_sd <- as.numeric(rep(x = NA, length(segments_gr)))
  # mcols(segments_gr)[unique(mcols(gccorrlogr_gr)$seg), "rho_logr_sd"] <- by(data = mcols(gccorrlogr_gr)$logr, INDICES = mcols(gccorrlogr_gr)$seg, FUN = sd, simplify = T, na.rm = T)
  # plot(mcols(segments_gr)$rho_logr, mcols(segments_gr)$rho_logr_sd)
  
  bafsegoverlaps <- findOverlaps(query = segments_gr, subject = phasedbaf)
  
  mcols(segments_gr)[unique(queryHits(bafsegoverlaps)), "nhetsnps"] <- by(data = subjectHits(bafsegoverlaps), INDICES = queryHits(bafsegoverlaps), FUN = length, simplify = T)
  mcols(segments_gr)[unique(queryHits(bafsegoverlaps)), "mu_baf"] <- by(data = mcols(phasedbaf)[subjectHits(bafsegoverlaps), "BAFphased"], INDICES = queryHits(bafsegoverlaps), FUN = mean, simplify = T, na.rm = T)
  mcols(segments_gr)[unique(queryHits(bafsegoverlaps)), "rho_baf"] <- by(data = mcols(phasedbaf)[subjectHits(bafsegoverlaps), "BAFphased"], INDICES = queryHits(bafsegoverlaps), FUN = sd, simplify = T, na.rm = T)
  
  mcols(segments_gr)$mu_vaf <- mcols(segments_gr)$mu_baf - (1-rhopsi$purity)/((2*(1-rhopsi$purity)+rhopsi$purity*rhopsi$ploidy)*2^mcols(segments_gr)$mu_logr)
  # unlikely to call variants present at a VAF of < 5%, so set to BAF so unlikely to call anything.
  mcols(segments_gr)$mu_vaf <- ifelse(mcols(segments_gr)$mu_vaf < .05, mcols(segments_gr)$mu_baf, mcols(segments_gr)$mu_vaf)
  
  write_tsv(x = as.data.frame(segments_gr), file = file.path(sampledir, paste0(sampleid, "_baf_logr_summarised_segments.txt")))
  
  return(segments_gr)
}


get_snv_mnvs <- function(sampleid, snvdir) {
  snv_mnvfile <- list.files(path = snvdir, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz$"), full.names = T, recursive = T)
  if (file.exists(snv_mnvfile)) {
    som_vcf <- readVcf(file = snv_mnvfile)
  } else {
    return(NULL)
  }
  # drop variants without allele counts reported
  #### MOD here to NOT BOTH NA and fix missing read counts by setting to 0
  som_vcf <- som_vcf[which(!(is.na(info(som_vcf)$t_alt_count) & is.na(info(som_vcf)$t_ref_count)))]
  som <- rowRanges(som_vcf)
  mcols(som) <- cbind(mcols(som), info(som_vcf))
  som$ALT <- unlist(som$ALT)
  som$t_alt_count[is.na(som$t_alt_count)] <- 0
  som$t_ref_count[is.na(som$t_ref_count)] <- 0
  # snvhits <- findOverlaps(query = segments_gr, subject = som_vcf)
  # som_vcf <- som_vcf[subjectHits(snvhits)] # not all are within BAF segments, male X has large part removed
  # mcols(som_vcf)$seg <- queryHits(snvhits)
  return(som)
}




get_consensus_cn <- function(sampleid, cndir, bsgenome) {
  breaksfile <- file.path(CNDIR, paste0(sampleid, ".consensus.20170119.somatic.cna.annotated.txt"))
  if (!file.exists(breaksfile))
    return(NULL)
  breaks <- GRanges(read.delim(file = breaksfile, as.is = T), seqinfo = seqinfo(bsgenome))
  # breaks_gr <- GRanges(seqnames = breaks$chromosome, ranges = IRanges(start = breaks$start, end = breaks$end), seqinfo = seqinfo(bsgenome))
  return(breaks)
}



test_clean_sites <- function(sampledir, sampleid, segments_gr, phasedbaf_gr, snvs_vcf, bsgenome, logr, ncores = 12, presched = T, pseudocountrange = c(50, 1000), subsample_optim = T, recompute_pseudocoverage = T, immune_loci = NA, germline_svs = NA) {

  ## drop segment if < 100 hetSNPs / vector
  # segments_gr <- segments_gr[mcols(segments_gr)$nhetsnps >= minsegmentSNPcount]
  
  ## drop SNVs and hetSNPs that overlap
  mutsnphits <- findOverlaps(query = phasedbaf_gr, subject = snvs_vcf)
  if (length(mutsnphits) > 0) {
    phasedbaf_gr <- phasedbaf_gr[-queryHits(mutsnphits)]
    snvs_vcf <- snvs_vcf[-subjectHits(mutsnphits)]
  }
  
  # get idxs of SNPs in (long) segments  
  phasedbaf_gr <- phasedbaf_gr[which(phasedbaf_gr %within% segments_gr)]
  
  # filter out immune regions
  phasedbaf_gr <- subsetByOverlaps(x = phasedbaf_gr, ranges = immune_loci, invert = T)
  
  snpseghits <- findOverlaps(query = phasedbaf_gr, subject = segments_gr)
  phasedbaf_gr <- phasedbaf_gr[queryHits(snpseghits)]
  
  # create vectors for repeated use below
  bbparams <- data.frame(size = mcols(phasedbaf_gr)$refCount + mcols(phasedbaf_gr)$altCount, 
                         shape1 = mcols(segments_gr)[subjectHits(snpseghits), "mu_baf"], 
                         shape2 = 1-mcols(segments_gr)[subjectHits(snpseghits), "mu_baf"],
                         majcount = mcols(phasedbaf_gr)$MajCount,
                         maxcount = pmax(mcols(phasedbaf_gr)$refCount, mcols(phasedbaf_gr)$altCount))
  
  # sampleidxs <- round(seq(1, nrow(lmdf), length.out = min(nrow(lmdf), 1e5)))
  
  pseudocoveragefile <- file.path(sampledir, paste0(sampleid, "_pseudocount_calibrated.txt"))
  
  if (!recompute_pseudocoverage && file.exists(pseudocoveragefile)) {
    pseudocount_calibr <- read.delim(file = pseudocoveragefile, as.is = T, header = F)[1,1]
  } else {
    if (subsample_optim) {
      subidxs <- seq(from = 1, to = nrow(bbparams), by = max(round( nrow(bbparams)/1e5 ), 1))
      pseudocount_calibr <- optimise(f = get_model_slope_dev, bbpar = bbparams[subidxs, ], ncores = ncores, presched = presched, interval = pseudocountrange)$minimum
    } else {
      pseudocount_calibr <- optimise(f = get_model_slope_dev, bbpar = bbparams, ncores = ncores, presched = presched, interval = pseudocountrange)$minimum
    }
  }
  
  # modelfits <- calibrate_model(pcounts = pseudocounts, bbpar = bbparams, ncores = ncores, presched = presched)
  
  ## QC of the fit below, not used further
  mcols(phasedbaf_gr)$pvalphased <- mcmapply(size = bbparams$size,
                     q = bbparams$majcount,
                     shape1 = bbparams$shape1*pseudocount_calibr,
                     shape2 = bbparams$shape2*pseudocount_calibr,
                     FUN = betabinom.test.ab,
                     MoreArgs = list(alternative = "two.sided"),
                     mc.preschedule = presched,
                     mc.cores = ncores)

  ggd.qqplot(sampleid = sampleid, sampledir = sampledir, suffix = "_phasedSNPs_QQplot", pvector = mcols(phasedbaf_gr)$pvalphased)
  
  # write phased + tested BAF to disk
  write.table(x = as.data.frame(phasedbaf_gr)[, c("seqnames", "start", "end", "refCount", "altCount", "MajCount", "pvalphased")], file = file.path(sampledir, paste0(sampleid, "_phasedbaf_tested.txt")), sep = "\t", quote = F, col.names = T, row.names = F)
  

  ######## somatic variants
  # tag immune regions
  mcols(snvs_vcf)$immune_locus <- snvs_vcf %within% immune_loci
  
  # subset SNVs to those on (long) segments (and PAR when appriopriate)
  snvseghits <- findOverlaps(query = snvs_vcf, subject = segments_gr)
  snvs_vcf <- snvs_vcf[queryHits(snvseghits)]
  
  if (!is.na(germline_svs)) {
    mcols(snvs_vcf)$germline_sv <- overlapsAny(query = snvs_vcf, subject = germline_svs)
  } else {
    mcols(snvs_vcf)$germline_sv <-  F
  }
  
  # independent filtering (check to see whether useful at PCAWG coverage)
    pfilter <- mcmapply(size = mcols(snvs_vcf)$t_alt_count + mcols(snvs_vcf)$t_ref_count,
                        q = mcols(snvs_vcf)$t_alt_count + mcols(snvs_vcf)$t_ref_count,
                        shape1 = mcols(segments_gr)[subjectHits(snvseghits), "mu_vaf"]*pseudocount_calibr,
                        shape2 = (1-mcols(segments_gr)[subjectHits(snvseghits), "mu_vaf"])*pseudocount_calibr,
                        FUN = betabinom.test.ab,
                        MoreArgs = list(alternative = "greater"),
                        mc.preschedule = presched,
                        mc.cores = ncores)
  ### end of filtering
  
  
  pvalsnv <- mcmapply(size = mcols(snvs_vcf)$t_alt_count + mcols(snvs_vcf)$t_ref_count,
                     q = mcols(snvs_vcf)$t_alt_count,
                     shape1 = mcols(segments_gr)[subjectHits(snvseghits), "mu_vaf"]*pseudocount_calibr,
                     shape2 = (1-mcols(segments_gr)[subjectHits(snvseghits), "mu_vaf"])*pseudocount_calibr,
                     FUN = betabinom.test.ab,
                     MoreArgs = list(alternative = "greater"),
                     mc.preschedule = presched,
                     mc.cores = ncores)
  
  pfilter[which(pfilter < .Machine$double.eps)] <- .Machine$double.eps
  pvalsnv[which(pvalsnv < .Machine$double.eps)] <- .Machine$double.eps
  # padjsnv <- p.adjust(pvalsnv, method = "fdr")
  
  ggd.qqplot(sampleid = sampleid, sampledir = sampledir, suffix = "_SNVs_QQplot", pvector = pvalsnv)
  
  # get slope for SNV QQ-plot
  lmdf_snv <- data.frame(exp = -log10(1:length(pvalsnv)/length(pvalsnv)), obs = -log10(sort(pvalsnv, decreasing = F, na.last = T)))
  lmdf_snv <- lmdf_snv[is.finite(lmdf_snv$obs),]
  snv_slope <- coef(MASS::rlm(formula = obs ~ 0 + exp, data = lmdf_snv))
  

  write.table(x = data.frame(pseudocount_calibr, snv_slope), file = pseudocoveragefile, sep = "\t", quote = F, col.names = F, row.names = F)
  
  # and associate with "corrected, somatic" BAF
  # info(header(snvs_vcf)) <- rbind(info(header(snvs_vcf)), DataFrame(Number = "1", Type = "Float", Description = "Corrected, somatic BAF of segment", row.names = "VAFseg"))
  # info(snvs_vcf)$VAFseg <- mcols(segments_gr)[subjectHits(snv_segment_hits), "mu_vaf"]
  # mcols(snvs_vcf)$VAFseg <- mcols(segments_gr)[subjectHits(snv_segment_hits), "mu_vaf"]
  # pvalsnv <- test_betabin_model(size = info(snvs_vcf)$t_alt_count + info(snvs_vcf)$t_ref_count,
  #                               q = info(snvs_vcf)$t_alt_count, 
  #                               shape1 = mcols(snvs_vcf)$VAFseg*1000, 
  #                               shape2 = (1 - mcols(snvs_vcf)$VAFseg)*1000,
  #                               idxs = 1:length(snvs_vcf), alternative = "greater")
  
  # QC
  # hist(pphased$pval)
  # plot(1:nrow(pphased), log10(pphased$pval))
  # plot(1:nrow(pphased), log10(pphased$padj))
  # ggd.qqplot(pvector = pphased$pval)
  # ggd.qqplot(pvector = pvalsnv)
  # ggd.qqplot(pvector = padjsnv)
  # ggd.qqplot(pvector = pphased$padj)
  # ggd.qqplot(pvector = pphased$pvalre)
  # ggd.qqplot(pvector = pphased$padjre)
  # 
  # qqplot2(pvector = pphased$pval, pvector_conserv = pphased$pvalre)
  
  # snppassidxs <- which(padjre <= fwer)
  # snvpassidxs <- which(padjsnv <= fwer)
  # 
  # phasedbaf_gr <- phasedbaf_gr[snppassidxs]
  # snvs_vcf <- snvs_vcf[snvpassidxs]
  
  mcols(snvs_vcf)$pval <- pvalsnv
  mcols(snvs_vcf)$pfilt <- pfilter
  mcols(snvs_vcf) <- cbind(mcols(snvs_vcf), mcols(segments_gr)[subjectHits(snvseghits), c("total_cn", "major_cn", "minor_cn") ])
  
  # add number of proximal hetSNPs (biases allele counts when 2x alt on ref/alt allele)
  mcols(snvs_vcf)$nhetsnps25bp <- countOverlaps(query = snvs_vcf, subject = phasedbaf_gr, maxgap = 25)
  
  # annotate hits with upstream and downstream BAF/LogR of SNPS
  # BAF
  preidxs <- precede(x = snvs_vcf, subject = phasedbaf_gr)
  postidxs <- follow(x = snvs_vcf, subject = phasedbaf_gr)
  mcols(snvs_vcf)$bafpos_pre <- start(phasedbaf_gr)[postidxs]
  mcols(snvs_vcf)$bafpos_post <- start(phasedbaf_gr)[preidxs]
  mcols(snvs_vcf)$bafpval_pre <- mcols(phasedbaf_gr)$pvalphased[postidxs]
  mcols(snvs_vcf)$bafpval_post <- mcols(phasedbaf_gr)$pvalphased[preidxs]
  
  # Log R
  logrsub <- logr[na.omit(unique(c(precede(x = snvs_vcf, subject = logr), follow(x = snvs_vcf, subject = logr))))]
  logrsegidxs <- nearest(x = logrsub, subject = segments_gr)

  mcols(logrsub)$pval <- pnorm(q = mcols(logrsub)$logr, mean = mcols(segments_gr)[logrsegidxs, "mu_logr"],
                                      sd = mcols(segments_gr)[logrsegidxs, "rho_logr"])
  mcols(logrsub)$pval <- 2*ifelse( mcols(logrsub)$pval >= .5, 1 - mcols(logrsub)$pval, mcols(logrsub)$pval)
  
  preidxs <- precede(x = snvs_vcf, subject = logrsub)
  postidxs <- follow(x = snvs_vcf, subject = logrsub)
  mcols(snvs_vcf)$logrpos_pre <- start(logrsub)[postidxs]
  mcols(snvs_vcf)$logrpos_post <- start(logrsub)[preidxs]
  mcols(snvs_vcf)$logrpval_pre <- mcols(logrsub)[postidxs, "pval"]
  mcols(snvs_vcf)$logrpval_post <- mcols(logrsub)[preidxs, "pval"]

  
  ## filter out missed segments: identify as multiple hits within certain window
  # allhits <- GRanges(seqnames = c(seqnames(hitsnps_gr), seqnames(hitsnvs_gr)), ranges = c(ranges(hitsnps_gr), ranges(rowRanges(hitsnvs_gr))) )
  # get those loci which may have hits within 1kb but no more > 1kb and <= 10 kb
  # gc_candidate_idxs <- which(countOverlaps(query = allhits, subject = allhits, maxgap = testwindow/2) - 
  #                              countOverlaps(query = allhits, subject = allhits, maxgap = conversionlength) == 0)
  # browser()

  outdf <- data.frame(chr = seqnames(snvs_vcf), start = start(snvs_vcf), end = end(snvs_vcf), ref = as.character(snvs_vcf$REF), alt = as.character(snvs_vcf$ALT))
  outcols <- c("VAF", "t_alt_count", "t_ref_count", "snv_near_indel", "Variant_Classification", "immune_locus", "germline_sv", "pval", "pfilt",
               "nhetsnps25bp", "total_cn", "major_cn", "minor_cn", "bafpos_pre", "bafpos_post", "bafpval_pre", "bafpval_post", 
               "logrpos_pre", "logrpos_post", "logrpval_pre", "logrpval_post",
               "1000genomes_AF", "n_ref_count", "n_alt_count", "n_total_cov",
               "Validation_status", "hg38clean")
  if (any(!outcols %in% colnames(mcols(snvs_vcf)))) {
    outcols <- intersect(x = outcols, y = colnames(mcols(snvs_vcf)))
    print("Subsetting output columns to omit hg38 checks etc")
  }
  outdf <- cbind(outdf, as.data.frame(mcols(snvs_vcf)[, outcols]))
  write.table(x = outdf, file = file.path(sampledir, paste0(sampleid, "_snv_mnv_infSites_annotated.txt")), sep = "\t", quote = F, col.names = T, row.names = F)
  
  return(list(hitvariants = snvs_vcf, pseudocount_calibr = pseudocount_calibr, snv_slope = snv_slope))
}


# calibrate_model <- function(pcounts, bbpar, ncores, presched) {
#   
#   lmdf <- data.frame(exp = -log10(1:nrow(bbpar)/nrow(bbpar)), obs = numeric(length = nrow(bbpar)))
#   
#   modelfits <- sapply(X = pcounts, FUN = get_model_slope, lmdf = lmdf, bbpar = bbpar, ncores = ncores, presched = presched, simplify = T)
#   
#   # take the model for which the coefficient's CI contains 1 and which is closest to 1
#   return(modelfits)
# }
# 


get_model_slope_dev <- function(pseudocount, bbpar, ncores, presched) {
  pvalphased <- mcmapply(size = bbpar$size,
                         q = bbpar$majcount,
                         shape1 = bbpar$shape1*pseudocount,
                         shape2 = bbpar$shape2*pseudocount,
                         FUN = betabinom.test.ab,
                         MoreArgs = list(alternative = "two.sided"),
                         mc.preschedule = presched,
                         mc.cores = ncores)
  
  lmdf <- data.frame(exp = -log10(1:nrow(bbpar)/nrow(bbpar)), obs = -log10(sort(pvalphased, decreasing = F)))
  lmdf <- lmdf[is.finite(rowSums(lmdf)),]
  qqlm <- MASS::rlm(formula = obs ~ 0 + exp, data = lmdf)
  
  # robci <- confint.default(object = qqlm, parm = "exp", level = 0.95)
  # outv <- c(pseudo = pseudocount, est = coef(qqlm)[[1]], lower = robci[1,1], upper = robci[1,2])
  slopedev <- (coef(qqlm) - 1)^2
  
  return(slopedev)
}



 annotate_gc_hit <- function(hitlocus, pseudocount_calibr, sampledir, sampleid, baf, logr, snvs, seg_gr, conversionlength = 1e3, testwindow = 2e4, fdr = .05, plotting = F) {
  hitseg <- seg_gr[which(seg_gr %over% hitlocus)]
  
  bafneighbours <- subsetByOverlaps(x = baf, ranges = hitlocus, maxgap = testwindow/2)
  # remove neighbours that sit on a different segment than the hit
  bafneighbours <- bafneighbours[bafneighbours %over% hitseg]
  
  hetsnpneighbours <- length(bafneighbours)
  # if (mcols(hitlocus)$type == "snp")
  #   hetsnpneighbours <- hetsnpneighbours - 1
  #   hitidx <- which(hitlocus %over% bafneighbours)
  # } else {
  #   
  # }
  
  if (hetsnpneighbours == 0) {
    mcols(hitlocus) <- data.frame(mcols(hitlocus),
                                  hetsnpneighbours = 0,
                                  mininfosnps = 0,
                                  nsnphits = 0,
                                  nsnphitsoutofrange = 0,
                                  predists = as.character(NA),
                                  postdists = as.character(NA),
                                  prepvals = as.character(NA),
                                  postpvals = as.character(NA),
                                  nsnvhits = 0,
                                  nsnvhitsoutofrange = 0,
                                  presnvdists = as.character(NA),
                                  postsnvdists = as.character(NA),
                                  presnvpvals = as.character(NA),
                                  postsnvpvals = as.character(NA),
                                  nlogrneighbours = 0,
                                  nlogrhits = 0,
                                  p_empirical = 0)
    return(hitlocus)
  }
  
  # if (length(bafneighbours) - 1 <= minothersnps)
  #   return(NULL)
  
  # minimal number of snps on either side of the hit? assures we have info on either side
  mininfosnps <- min(sum(start(bafneighbours) < start(hitlocus)), sum(start(bafneighbours) > start(hitlocus)))
  
  # hitidx <- nearest(x = hitlocus, subject = bafneighbours, select = "arbitrary")
  # is_terminal <- hitidx == 1 || hitidx == length(bafneighbours)
  # if (hitidx == 1 || hitidx == length(bafneighbours))
  #   return(NULL)
  
  # test BAF of neighbours
  # pphased_sub <- test_betabin_model(size = mcols(bafneighbours)$refCount + mcols(bafneighbours)$altCount,
  #                                   q = mcols(bafneighbours)$MajCount, 
  #                                   shape1 = mcols(bafneighbours)$BAFseg*1000, 
  #                                   shape2 = (1 - mcols(bafneighbours)$BAFseg)*1000,
  #                                   idxs = 1:length(bafneighbours))
  
  pvalphased <- mapply(size = mcols(bafneighbours)$refCount + mcols(bafneighbours)$altCount,
                     q = mcols(bafneighbours)$MajCount,
                     FUN = betabinom.test.ab,
                     MoreArgs = list(alternative = "two.sided",
                                     shape1 = mcols(hitseg)$mu_baf*pseudocount_calibr,
                                     shape2 = (1-mcols(hitseg)$mu_baf)*pseudocount_calibr))
  
  
  padjphased <- p.adjust(pvalphased, method = "fdr")
  
  nsnphits <- sum(padjphased <= fdr)
  if (mcols(hitlocus)$type == "snp")
    nsnphits <- nsnphits - 1
  nsnphitsoutofrange <- sum(!overlapsAny(query = bafneighbours[padjphased <= fdr], subject = hitlocus, maxgap = conversionlength))
  
  snpdists <- start(bafneighbours) - start(hitlocus)
  predists <- paste0(rev(snpdists[snpdists < 0]), collapse = ";")
  postdists <- paste0(snpdists[snpdists > 0], collapse = ";")
  
  prepvals <- paste0(rev(pvalphased[snpdists < 0]), collapse = ";")
  postpvals <- paste0(pvalphased[snpdists > 0], collapse = ";")
  
  
  #### snvs
  snvneighbours <- subsetByOverlaps(x = snvs, ranges = hitlocus, maxgap = testwindow/2)
  # remove neighbours that sit on a different segment than the hit
  snvneighbours <- snvneighbours[snvneighbours %over% hitseg]
  
  nsnvneighbours <- length(snvneighbours)
  if (mcols(hitlocus)$type == "snv") {
    nsnvneighbours <- nsnvneighbours - 1
  }
  
  if (nsnvneighbours > 0) {
    pvalsnv <- mapply(size = info(snvneighbours)$t_alt_count + info(snvneighbours)$t_ref_count,
                         q = info(snvneighbours)$t_alt_count,
                         FUN = betabinom.test.ab,
                         MoreArgs = list(alternative = "greater",
                                         shape1 = mcols(hitseg)$mu_vaf*pseudocount_calibr,
                                         shape2 = (1-mcols(hitseg)$mu_vaf)*pseudocount_calibr))
    
    pvalsnv[pvalsnv < .Machine$double.eps] <- .Machine$double.eps
    
    padjsnv <- p.adjust(pvalsnv, method = "fdr")
    
    nsnvhits <- sum(padjsnv <= fdr)
    if (mcols(hitlocus)$type == "snv")
      nsnvhits <- nsnvhits - 1
    nsnvhitsoutofrange <- sum(!overlapsAny(query = snvneighbours[padjsnv <= fdr], subject = hitlocus, maxgap = conversionlength))
    
    snvdists <- start(snvneighbours) - start(hitlocus)
    presnvdists <- paste0(rev(snvdists[snvdists < 0]), collapse = ";")
    postsnvdists <- paste0(snvdists[snvdists > 0], collapse = ";")
    
    presnvpvals <- paste0(rev(pvalsnv[snvdists < 0]), collapse = ";")
    postsnvpvals <- paste0(pvalsnv[snvdists > 0], collapse = ";")
  } else {
    nsnvhits <- 0
    nsnvhitsoutofrange <- 0
    presnvdists <- ""
    postsnvdists <- ""
    presnvpvals <- ""
    postsnvpvals <- ""
  }
  
  
  # OK if only hits are within 1kb of SNP
  # all_hits_in_range <- all(which(pphased_sub$padj <= fdr) %in% subjectHits(findOverlaps(query = hitlocus, subject = bafneighbours, maxgap = conversionlength)))
  # if (!all_hits_in_range)
  #   return(NULL)
  
  # check whether the rephased model is considerably better than the original, or whether phasing was correct.
  # prephased <- test_betabin_model(size = mcols(hitlocus)$refCount + mcols(hitlocus)$altCount,
  #                                   q = max(mcols(hitlocus)$refCount, mcols(hitlocus)$altCount), 
  #                                   shape1 = mcols(hitlocus)$BAFseg*1000, 
  #                                   shape2 = (1 - mcols(hitlocus)$BAFseg)*1000,
  #                                   idxs = 1)
  
  
  
  ## incorporate similar LogR checks
  # get all loci within 10kb either side
  logrneighbours <- subsetByOverlaps(x = logr, ranges = hitlocus, maxgap = testwindow/2)
  # remove neighbours that sit on a different segment than the hit
  logrneighbours <- logrneighbours[logrneighbours %over% hitseg]
  nlogrneighbours <- length(logrneighbours)
  
  mcols(logrneighbours)$pval <- pnorm(q = mcols(logrneighbours)$logr, mean = mcols(hitseg)$mu_logr,
                                      sd = mcols(hitseg)$rho_logr)
  mcols(logrneighbours)$pval <- 2*ifelse( mcols(logrneighbours)$pval >= .5,
                                          1 - mcols(logrneighbours)$pval,
                                          mcols(logrneighbours)$pval)
  
  # check how many more extreme than hit
  if (mcols(hitlocus)$type == "snp") {
    p_empirical <- (sum(mcols(subsetByOverlaps(x = logrneighbours, ranges = hitlocus))$pval >= mcols(logrneighbours)$pval) - 1) / (length(logrneighbours) - 1)
  } else {
    # check closest SNP if hit is SNV
    p_empirical <- (sum(mcols(logrneighbours[nearest(x = hitlocus, subject = logrneighbours)])$pval >= mcols(logrneighbours)$pval) - 1) / (length(logrneighbours) - 1)
  }
  # if (p_empirical <= .05)
  #   return(NULL)
  
  # test if all derive from normal distribution 
  # df = 2*var/(var-1) (for t-distrib, )
  
  # logrpvals <- data.frame(chr = seqnames(logrneighbours), pos = start(logrneighbours), logr = mcols(logrneighbours)$logr, pval = logrpvals)
  
  # from local loci, take maximally informative set + all within 1kb
  logrintervals <- GRanges(seqnames = seqnames(hitlocus), IRanges(start = seq(from = start(hitlocus)-testwindow/2, to =start(hitlocus)+testwindow/2, by = 1000), width = 1), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
  logrneighbours_spread <- unique(c(logrneighbours[nearest(x = logrintervals, subject = logrneighbours, select = "arbitrary")],
                                    subsetByOverlaps(x = logrneighbours, ranges = hitlocus, maxgap = conversionlength/2)))
  
  mcols(logrneighbours_spread)$padj <- p.adjust(mcols(logrneighbours_spread)$pval, method = "fdr")
  
  # formal testing OK only if no hits
  nlogrhits <- sum(mcols(logrneighbours_spread)$padj <= fdr)
  # if (logr_hits)
  #   return(NULL)
  
  mcols(hitlocus) <- data.frame(mcols(hitlocus),
                                hetsnpneighbours = hetsnpneighbours,
                                mininfosnps = mininfosnps,
                                nsnphits = nsnphits,
                                nsnphitsoutofrange = nsnphitsoutofrange,
                                predists = predists,
                                postdists = postdists,
                                prepvals = prepvals,
                                postpvals = postpvals,
                                nsnvhits = nsnvhits,
                                nsnvhitsoutofrange = nsnvhitsoutofrange,
                                presnvdists = presnvdists,
                                postsnvdists = postsnvdists,
                                presnvpvals = presnvpvals,
                                postsnvpvals = postsnvpvals,
                                nlogrneighbours = nlogrneighbours,
                                nlogrhits = nlogrhits,
                                p_empirical = p_empirical)
  
  is_clean_gc <- (hetsnpneighbours >= 3 & mininfosnps >= 1 & nsnphitsoutofrange <= 1 &
                    (as.integer(strsplit(x = predists, split = ";")[[1]][1]) >= -2000 | as.integer(strsplit(x = postdists, split = ";")[[1]][1]) <= 2000) &
                    nsnvhitsoutofrange == 0 & nlogrhits == 0 & p_empirical > .05)
  
  if (plotting & is_clean_gc) {
    # if all tests passed, plot
    plotdf <- data.frame(chr = seqnames(bafneighbours), pos = start(bafneighbours), VAF = mcols(bafneighbours)$BAFphased,
                         nsuccess = mcols(bafneighbours)$MajCount, nfail = mcols(bafneighbours)$MinCount)
    plotdf[, c("lower", "upper")] <- t(mapply(FUN = HDIofICDF,
                                              shape1 = plotdf$nsuccess + 1,
                                              shape2 = plotdf$nfail + 1,
                                              MoreArgs = list(ICDFname = qbeta, credMass = .95)))
    plotdf$type <- "snp"
    # pphased_sub$VAF <- mcols(bafneighbours)$BAFphased
    # pphased_sub$pos <- start(bafneighbours)
    
    if (length(snvneighbours) > 0) {
      plotdfsnv <- data.frame(chr = seqnames(snvneighbours), pos = start(snvneighbours), VAF = info(snvneighbours)$VAF,
                              nsuccess = info(snvneighbours)$t_alt_count, nfail = info(snvneighbours)$t_ref_count)
      plotdfsnv[, c("lower", "upper")] <- t(mapply(FUN = HDIofICDF,
                                                   shape1 = plotdfsnv$nsuccess + 1,
                                                   shape2 = plotdfsnv$nfail + 1,
                                                   MoreArgs = list(ICDFname = qbeta, credMass = .95)))
      plotdfsnv$type <- "snv"
      plotdf <- rbind(plotdf, plotdfsnv)
    }
    plotdf$type <- factor(x = plotdf$type, levels = c("snp", "snv"))
    
    hlinepars <- c(mcols(hitseg)$mu_baf, HDIofICDF(ICDFname = qbeta, credMass = .95, shape1 = mcols(hitseg)$mu_baf*pseudocount_calibr, shape2 = (1 - mcols(hitseg)$mu_baf)*pseudocount_calibr))
    hline2pars <- c(mcols(hitseg)$mu_vaf, HDIofICDF(ICDFname = qbeta, credMass = .95, shape1 = mcols(hitseg)$mu_vaf*pseudocount_calibr, shape2 = (1 - mcols(hitseg)$mu_vaf)*pseudocount_calibr))
    
    p1 <- ggplot(data = plotdf, mapping = aes(x = pos, y = VAF)) + coord_cartesian(ylim = c(0,1))
    p1 <- p1 + geom_pointrange(mapping = aes(ymin = lower, ymax = upper, colour = type), alpha = .7, show.legend = F)
    p1 <- p1 + geom_hline(yintercept = hlinepars[1]) + geom_ribbon(mapping = aes(x = seq(from = start(hitlocus)-10000, to = start(hitlocus)+10000, length.out = nrow(plotdf)), ymin = hlinepars[2], ymax = hlinepars[3]), fill = "red", alpha = .2)
    p1 <- p1 + geom_hline(yintercept = hline2pars[1]) + geom_ribbon(mapping = aes(x = seq(from = start(hitlocus)-10000, to = start(hitlocus)+10000, length.out = nrow(plotdf)), ymin = hline2pars[2], ymax = hline2pars[3]), fill = "blue", alpha = .2)
    p1 <- p1 + theme_minimal() + labs(title = paste0("chr", seqnames(hitlocus), ":", start(hitlocus)), x = "position", y = "phased BAF") + theme(plot.title = element_text(hjust=0.5))
    # p1

    p2 <- ggplot(data = as.data.frame(logrneighbours), mapping = aes(x = start, y = logr)) + geom_point()
    p2 <- p2 + geom_vline(xintercept = start(hitlocus))
    p2 <- p2 + geom_hline(yintercept = mcols(hitseg)$mu_logr)
    p2 <- p2 + geom_ribbon(mapping = aes(x = seq(from = start(hitlocus)-testwindow/2, to = start(hitlocus)+testwindow/2, length.out = length(logrneighbours)),
                                         ymin = mcols(hitseg)$mu_logr - 2*mcols(hitseg)$rho_logr, 
                                         ymax = mcols(hitseg)$mu_logr + 2*mcols(hitseg)$rho_logr), fill = "grey50", alpha = .2)
    p2 <- p2 + theme_minimal() + labs(title = paste0("chr", seqnames(hitlocus), ":", start(hitlocus)), x = "position", y = "logR") + theme(plot.title = element_text(hjust=0.5))
    p2 <- p2 + coord_cartesian(ylim = c(-2,2))
    # p2
    
    ggsave(filename = file.path(sampledir, "figures", paste0(sampleid, "_chr", seqnames(hitlocus), "_", start(hitlocus), "_bafplot.png")), plot = p1)
    ggsave(filename = file.path(sampledir, "figures", paste0(sampleid, "_chr", seqnames(hitlocus), "_", start(hitlocus), "_logrplot.png")), plot = p2)
    }
  
  return(hitlocus)
}


annotate_gc_hits <- function(hitvariants, pseudocount_calibr, sampledir, sampleid, baf, logr, snvs, seg_gr, conversionlength = 1e3, testwindow = 2e4, fdr = .05, plotting = T) {
  
  if (length(hitvariants) == 0 | is.null(hitvariants)) {
    write_tsv(x = data.frame(seqnames = character(),
                             start = integer(),
                             end = integer(),
                             width = integer(),
                             strand = character(),
                             type = character(),
                             pval = numeric(),
                             padj = numeric(),
                             hetsnpneighbours = integer(),
                             mininfosnps = integer(),
                             nsnphits = integer(),
                             nsnphitsoutofrange = integer(),
                             predists = character(),
                             postdists = character(),
                             prepvals = character(),
                             postpvals = character(),
                             nsnvhits = integer(),
                             nsnvhitsoutofrange = integer(),
                             presnvdists = character(),
                             postsnvdists = character(),
                             presnvpvals = character(),
                             postsnvpvals = character(),
                             nlogrneighbours = integer(),
                             nlogrhits = integer(),
                             p_empirical = numeric()), file = file.path(sampledir, paste0(sampleid, "_baflogr_gc_results.txt")))
    print(paste0("Sample ", sampleid, ": no conversion events detected"))
    return(hitvariants)
  }
  
  if (plotting)
    dir.create(file.path(sampledir, "figures"), showWarnings = F)
  
  hitvariants_grlist_annotated <- lapply(X = split(hitvariants, f = 1:length(hitvariants)), FUN = annotate_gc_hit, pseudocount_calibr = pseudocount_calibr, sampledir = sampledir, sampleid = sampleid, baf = baf, logr = logr, seg_gr = seg_gr, snvs = snvs, plotting = plotting, conversionlength = conversionlength, testwindow = testwindow, fdr = fdr)
  hitvariants_gr_annotated <- unlist(GRangesList(hitvariants_grlist_annotated))
  
  write_tsv(x = as.data.frame(hitvariants_gr_annotated), file = file.path(sampledir, paste0(sampleid, "_baflogr_gc_results.txt")))
  
  # hitsnvs_gr_annotated <- hitvariants_gr_annotated[mcols(hitvariants_gr_annotated)$type == "snv"]
  snvmatches <- findOverlaps(query = hitvariants_gr_annotated, subject = snvs)
  
  # snvs <- snvs[subjectHits(snvmatches)]
  addheader <- DataFrame(Number = rep("1", 11), Type = rep("Float", 11), Description = paste0("Gene conversion annotation ", 1:11), row.names = colnames(mcols(hitvariants_gr_annotated))[c( 2:7, 12:13, 18:20)])
  info(header(snvs)) <- rbind(info(header(snvs)), addheader)
  
  info(snvs)$pval <- as.numeric(NA)
  info(snvs)$padj <- as.numeric(NA)
  info(snvs)$hetsnpneighbours <- as.integer(NA)
  info(snvs)$mininfosnps <- as.integer(NA)
  info(snvs)$nsnphits <- as.integer(NA)
  info(snvs)$nsnphitsoutofrange <- as.integer(NA)
  info(snvs)$nsnvhits <- as.integer(NA)
  info(snvs)$nsnvhitsoutofrange <- as.integer(NA)
  info(snvs)$nlogrneighbours <- as.integer(NA)
  info(snvs)$nlogrhits <- as.integer(NA)
  info(snvs)$p_empirical <- as.numeric(NA)
  
  info(snvs)[subjectHits(snvmatches), colnames(mcols(hitvariants_gr_annotated))[c( 2:7, 12:13, 18:20)] ] <- mcols(hitvariants_gr_annotated)[queryHits(snvmatches), c( 2:7, 12:13, 18:20)]
  
  writeVcf(obj = snvs, filename = file.path(sampledir, paste0(sampleid, "_snvs_baflogr_gc_results_fullannot.vcf")))
  
  return(hitvariants_gr_annotated)
}


trim_XY_to_PAR <- function(segments_gr, bsgenome) {
  PAR <- GRanges(seqnames = c("X", "X", "Y", "Y"), ranges = IRanges(start = c(60001, 154931044, 10001, 59034050), end = c(2699520, 155270559, 2649520, 59373565)), seqinfo = seqinfo(bsgenome))
  PAR_segmented <- subsetByOverlaps(x = disjoin(c(segments_gr[seqnames(segments_gr) %in% c("X", "Y")], PAR)), ranges = PAR)
  segments_gr <- c(segments_gr[!seqnames(segments_gr) %in% c("X", "Y")], PAR_segmented)
  return(segments_gr)
}




call_parallel_violations <- function(sampleid, sampledir, phasingdir, nboot = 1000, alpha = .1) {
  
  # sampledir <- file.path(vafhitsdir, sampleid)
  
  # read data if exists, otherwise exit
  vafhitsfile <- file.path(sampledir, paste0(sampleid, "_snv_mnv_infSites_annotated.txt"))
  phasinghitsfile <- file.path(phasingdir, sampleid, paste0(sampleid, "_tumour_snv-snp_phased.txt.gz"))
  testedbaffile <- file.path(sampledir, paste0(sampleid, "_phasedbaf_tested.txt"))
  
  if (any(!file.exists(vafhitsfile, testedbaffile))) return(NULL)
  
  vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
  testedbaf <- read_tsv(file = testedbaffile, col_types = "ciiiiin")
  
  vafhitsidxs <- paste0(vafhitsdf$chr, "_", vafhitsdf$start)
  testedbaf$snpidx <- paste0(testedbaf$seqnames, "_", testedbaf$start)
  
  
  if (file.exists(phasinghitsfile)) {
    phasinghits <- read.delim(file = phasinghitsfile, as.is = T)
    # SNV-SNP phasing only directly informative
    phasinghits <- phasinghits[which(phasinghits$type1 != phasinghits$type2), ]
    
    # make sure there is a decent amount of pair coverage necessary to see anything (i.e. 2 reads for each SNP allele, 4 reads from the somatic variant allele, fair since we won't pick up anything subclonal anyhow)
    coveredphasinghits <- which(ifelse(phasinghits$type1 == "SNP",
                                       phasinghits$Num_ref_ref + phasinghits$Num_ref_alt >= 2 & phasinghits$Num_alt_ref + phasinghits$Num_alt_alt >= 2 & phasinghits$Num_ref_alt + phasinghits$Num_alt_alt >= 4 ,
                                       phasinghits$Num_ref_ref + phasinghits$Num_alt_ref >= 2 & phasinghits$Num_alt_alt + phasinghits$Num_ref_alt >= 2 & phasinghits$Num_alt_ref + phasinghits$Num_alt_alt >= 4))
    phasinghits <- phasinghits[coveredphasinghits, ]
    
    # create snv/snp indices for quick overlaps
    phasinghits$snvidx <- ifelse(phasinghits$type1 == "SNV", paste0(phasinghits$chr, "_", phasinghits$pos1),  paste0(phasinghits$chr, "_", phasinghits$pos2))
    phasinghits$snpidx <- ifelse(phasinghits$type1 == "SNV", paste0(phasinghits$chr, "_", phasinghits$pos2),  paste0(phasinghits$chr, "_", phasinghits$pos1))
    
    # subset to the ones which have been evaluated for ISA violations
    phasinghitssub <- phasinghits[which(phasinghits$snvidx %in% vafhitsidxs), ]
    # testedbaf <- testedbaf[which(testedbaf$snpidx %in% phasinghitssub$snpidx), ]
    # and call violations from phasing data
    phasinghitssub$totalpwcov <- rowSums(phasinghitssub[, c("Num_ref_ref", "Num_ref_alt", "Num_alt_alt", "Num_alt_ref")])
    phasinghitssub$isaviol <- ifelse(phasinghitssub$type1 == "SNP",
                                     phasinghitssub$Num_ref_alt > 1 & phasinghitssub$Num_alt_alt > 1 & phasinghitssub$Num_ref_alt/phasinghitssub$totalpwcov > .1 & phasinghitssub$Num_alt_alt/phasinghitssub$totalpwcov > .1,
                                     phasinghitssub$Num_alt_ref > 1 & phasinghitssub$Num_alt_alt > 1 & phasinghitssub$Num_alt_ref/phasinghitssub$totalpwcov > .1 & phasinghitssub$Num_alt_alt/phasinghitssub$totalpwcov > .1)
    phasinghitssub$phasedsnpval <- testedbaf$pvalphased[match(x = phasinghitssub$snpidx, table = testedbaf$snpidx)]

  } else {
    phasinghitssub <- data.frame(chr = character(), pos1 = integer(), ref1 = character(), alt1 = character(), type1 = character(),
                                 pos2 = integer(), ref2 = character(), alt2 = character(), type2 = character(),
                                 Num_ref_ref = integer(), Num_alt_alt = integer(), Num_alt_ref = integer(), Num_ref_alt = integer(),
                                 phasing = character(), snvidx = character(), snpidx = character(), totalpwcov = integer(), isaviol = logical(), phasedsnpval = numeric())
  }
  
  # assign calls
  vafhitsdf$is_phaseable <- vafhitsidxs %in% phasinghitssub$snvidx
  vafhitsdf$is_confirmed <- vafhitsidxs %in% names(which(c(by(data = phasinghitssub, INDICES = phasinghitssub$snvidx, FUN = function(x) any(x$isaviol & x$phasedsnpval > 1e-3)))))
  
  # vafhitsdf$phasedsnpval <- NA
  # vafhitsdf[vafhitsdf$is_phaseable, "phasedsnpval"] <- c(by(data = phasinghitssub$phasedsnpval, INDICES = phasinghitssub$snvidx, FUN = min))[vafhitsidxs[vafhitsdf$is_phaseable]]
  # vafhitsdf$phasedsnpval <- c(by(data = phasinghitssub$phasedsnpval, INDICES = phasinghitssub$snvidx, FUN = min))[vafhitsidxs]
  
  # plotting during devel
  # p1 <- ggplot(data = vafhitsdfsub, mapping = aes(y = pval, x = is_confirmed)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + scale_y_log10()
  # p1
  # p2 <- ggplot(data = vafhitsdfsub, mapping = aes(y = pval, x = is_confirmed)) + geom_jitter(mapping = aes(colour = nhetsnps25bp >= 2)) + scale_y_log10()
  # p2
  # p2 <- ggplot(data = vafhitsdfsub[which(vafhitsdfsub$pval < 1e-4), ], mapping = aes(y = pmin(bafpval_pre, bafpval_post), x = is_confirmed)) + geom_jitter(mapping = aes(colour = pmin(bafpval_pre, bafpval_post) > 1e-3)) + scale_y_log10()
  # p2
  # p3 <- ggplot(data = vafhitsdfsub[which(vafhitsdfsub$pval < 1e-4), ], mapping = aes(y = pmin(logrpval_pre, logrpval_post), x = is_confirmed)) + geom_jitter(mapping = aes(colour = pmin(logrpval_pre, logrpval_post) > 1e-3)) + scale_y_log10()
  # p3
  # p3 <- ggplot(data = vafhitsdfsub[which(vafhitsdfsub$pval < 1e-4), ], mapping = aes(y = pval, x = is_confirmed, colour = is.na(total_cn) | minor_cn == 0)) + geom_jitter() + scale_y_log10()
  # p3
  # p4 <- ggplot(data = vafhitsdfsub, mapping = aes(y = pval, x = is_confirmed)) + geom_boxplot(outlier.shape = NA) + geom_jitter(mapping = aes(colour = pmin(bafpval_pre, bafpval_post) > 1e-3 & pmin(logrpval_pre, logrpval_post) > 1e-3)) + scale_y_log10()
  # p4
  
  # filter out variants where adjacent (het)SNPs do not behave (BAF/LogR) according to segment
  # filter out variants in the immune regions (as allele frequencies may be messed up)
  # also filter out sites which have ≥ 2 hetSNPs within 25bp window as they considerably bias the allelecounts when phased
  vafhitsdf$bafpval_comb <- rep(1, nrow(vafhitsdf))
  goodbafidxs <- which(vafhitsdf$bafpval_pre > .Machine$double.eps & vafhitsdf$bafpval_post > .Machine$double.eps)
  vafhitsdf$bafpval_comb[goodbafidxs] <- apply(X = vafhitsdf[goodbafidxs, c("bafpval_pre", "bafpval_post")], MARGIN = 1, FUN = function(x) sumlog(p = x)$p)
  
  vafhitsdf$logrpval_comb <- rep(1, nrow(vafhitsdf))
  goodlogridxs <- which(vafhitsdf$logrpval_pre > .Machine$double.eps & vafhitsdf$logrpval_post > .Machine$double.eps)
  vafhitsdf$logrpval_comb[goodlogridxs] <- apply(X = vafhitsdf[goodlogridxs, c("logrpval_pre", "logrpval_post")], MARGIN = 1, FUN = function(x) sumlog(p = x)$p)
  
  vafhitsdf_clean <- vafhitsdf[which(!is.na(vafhitsdf$pval) & vafhitsdf$minor_cn > 0 &
                                       pmin(vafhitsdf$bafpval_pre, vafhitsdf$bafpval_post) > 1e-3 &
                                       vafhitsdf$bafpval_comb > 1e-2 &
                                       pmin(vafhitsdf$logrpval_pre, vafhitsdf$logrpval_post) > 1e-3 &
                                       vafhitsdf$logrpval_comb > 1e-2 &
                                       !vafhitsdf$germline_sv &
                                       vafhitsdf$pfilt < 1e-3 &
                                       # !vafhitsdf$immune_locus &
                                       vafhitsdf$nhetsnps25bp < 2), ]

  #summary stats
  # par_phasing_conf <- sum(vafhitsroc$is_confirmed)
  # browser()
  # note that only in 1+1 regions do we expect phasing and VAF pipeline to produce the exact same matches (except for low CCF subclonal second hits), still, push all in.
  perfmets_all <- get_metrics(vafhitsdf = vafhitsdf_clean, sampleid = sampleid, sampledir = sampledir, nboot = nboot, plotting = T, alpha = alpha)
  # perfmets_hetero <- get_metrics(vafhitsdf = vafhitsdf_clean, rocidxs = which(vafhitsdf_clean$is_phaseable & vafhitsdf_clean$minor_cn > 0), sampleid = sampleid, sampledir = sampledir, nboot = nboot)
  # perfmets_diploid <- get_metrics(vafhitsdf = vafhitsdf_clean, rocidxs = which(vafhitsdf_clean$is_phaseable & vafhitsdf_clean$major_cn == 1 & vafhitsdf_clean$minor_cn == 1), sampleid = sampleid, sampledir = sampledir, nboot = nboot)
  # perfmets_hetero <- setNames(object = perfmets_hetero, nm = paste(names(perfmets_all), sep = "_", "hetero"))
  # perfmets_diploid <- setNames(object = perfmets_diploid, nm = paste(names(perfmets_all), sep = "_", "dipl"))
  
  finalhits <- vafhitsdf_clean[which(vafhitsdf_clean$pval <= perfmets_all[["cutoff"]] | vafhitsdf_clean$is_confirmed), ]
  estimtotal <- setNames(object = quantile(x = rbetabinom.ab(n = 1e4, size = nrow(vafhitsdf_clean), shape1 = sum(vafhitsdf_clean$is_confirmed)+.001, 
                                                             shape2 = sum(vafhitsdf_clean$is_phaseable & !vafhitsdf_clean$is_confirmed)+.001), probs = c(.025,.5,.975)),
                         nm = c("lower", "med", "upper"))
  estimtotal_diploid <- setNames(object = quantile(x = rbetabinom.ab(n = 1e4, size = sum(vafhitsdf_clean$major_cn == 1 & vafhitsdf_clean$minor_cn == 1, na.rm = T),
                                                                     shape1 = sum(vafhitsdf_clean$major_cn == 1 & vafhitsdf_clean$minor_cn == 1 & vafhitsdf_clean$is_confirmed, na.rm = T)+.001,
                                                                     shape2 = sum(vafhitsdf_clean$major_cn == 1 & vafhitsdf_clean$minor_cn == 1 & vafhitsdf_clean$is_phaseable & !vafhitsdf_clean$is_confirmed, na.rm = T)+.001), probs = c(.025,.5,.975)),
                                 nm = c("lower_diploid", "med_diploid", "upper_diploid"))
  sumstats <- c(tot_testable = nrow(vafhitsdf_clean), tot_hetero = sum(vafhitsdf_clean$minor_cn > 0, na.rm = T), 
                tot_diploid = sum(vafhitsdf_clean$major_cn == 1 & vafhitsdf_clean$minor_cn == 1, na.rm = T), nparallel = nrow(finalhits),
                nparallel_hetero = sum(finalhits$minor_cn > 0, na.rm = T), nparallel_dipl = sum(finalhits$major_cn == 1 & finalhits$minor_cn == 1, na.rm = T),
                tot_phaseable = sum(vafhitsdf_clean$is_phaseable), npar_phased = sum(finalhits$is_confirmed), npar_vaf = sum(finalhits$pval <= perfmets_all[["cutoff"]], na.rm = T),
                estimtotal, estimtotal_diploid, perfmets_all)
  # sumstats <- c(tot_testable = nrow(vafhitsdf_clean), tot_hetero = sum(vafhitsdf_clean$minor_cn > 0, na.rm = T),
  #               tot_diploid = sum(vafhitsdf_clean$major_cn == 1 & vafhitsdf_clean$minor_cn == 1, na.rm = T), 
  #               nparallel = nrow(finalhits), nparallel_diploid = sum(finalhits$major_cn == 1 & finalhits$minor_cn == 1, na.rm = T),
  #               tot_phaseable = nrow(vafhitsdf_clean$is_phaseable), npar_phased = sum(vafhitsdf_clean$is_confirmed),
  #               estimtotal, estimtotal_diploid, perfmets_all, perfmets_hetero, perfmets_diploid)
  
  # writing output
  write.table(x = finalhits, file = file.path(sampledir, paste0(sampleid, "_snv_mnv_infSites_finalhits.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(x = vafhitsdf, file = file.path(sampledir, paste0(sampleid, "_snv_mnv_infSites+phasing_annotated.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(x = phasinghitssub, file = file.path(sampledir, paste0(sampleid, "_tumour_snv-snp_phased_InfSitesInformative.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(x = sumstats, file = file.path(sampledir, paste0(sampleid, "_InfSites_VAFpipeline_summarystats.txt")), sep = "\t", quote = F, col.names = F, row.names = T)
  
  return(finalhits)
}


get_metrics <- function(vafhitsdf_clean, sampleid, sampledir, nboot, alpha = .1, plotting = F) {
  
  if (nrow(vafhitsdf_clean) > 0 ){
    # fdrcutoff <- sum(p.adjust(vafhitsdf_clean$pval, method = "fdr") <= alpha , na.rm = T)*alpha/sum(!is.na(vafhitsdf_clean$pval))
    # phaseableidxs <- which(vafhitsdf_clean$is_phaseable)
    
    # if there are phasing-confirmed parallel hits
    if (plotting) {
      invisible(get_performance_metrics(df = vafhitsdf_clean, alpha = alpha, plotting = plotting, sampleid = sampleid, sampledir = sampledir))
    }
    bootout <- boot::boot(data = vafhitsdf_clean, statistic = function(x, i) {get_performance_metrics(df = x[i,], alpha = alpha, sampleid = sampleid, sampledir = sampledir, plotting = F)}, R = nboot)
    perfmets <- setNames(object = colMeans(bootout$t, na.rm = T), names(bootout$t0))
    perfmets[["cutoff"]] <- 10^-perfmets[["cutoff"]]
  } else {
    perfmets <- c(cutoff = NA, prec = NA, rec = NA)
  }
  
  return(perfmets)
}


get_pval_cutoff <- function(df, alpha = .1) {
  cutoff <- -log10(sum(p.adjust(df$pval, method = "fdr") <= alpha, 1 , na.rm = T)*alpha/sum(!is.na(df$pval), 1))
  return(cutoff)
}


get_performance_metrics <- function(df, alpha, plotting = F, sampleid, sampledir) {
  
  cutoff <- -log10(sum(p.adjust(df$pval, method = "fdr") <= alpha, 1 , na.rm = T)*alpha/sum(!is.na(df$pval), 1))
  phaseableidxs <- which(df$is_phaseable)
  
  if (sum(df[phaseableidxs, "is_confirmed"]) == 0) {
    optimperf <- c(cutoff = cutoff, prec = NA, rec = NA)
    return(optimperf)
  }
  
  infsitespred <- prediction(predictions = -log10(df[phaseableidxs, "pval"]), labels = df[phaseableidxs, "is_confirmed"])
  infsitesperf <- performance(prediction.obj = infsitespred, measure = "prec", x.measure = "rec")
  
  fdridx <- tail(which(infsitesperf@alpha.values[[1]] >= cutoff), n = 1)
  if (fdridx == 1) {
    prec <- infsitesperf@y.values[[1]][fdridx+1]
  } else {
    prec <- infsitesperf@y.values[[1]][fdridx]
  }
  rec <- infsitesperf@x.values[[1]][fdridx]
  
  optimperf <- c(cutoff = cutoff, prec = prec, rec = rec)
  
  if (plotting) {
    p6 <- ggplot(data = as.data.frame(t(optimperf)), mapping = aes(y = prec, x = rec))
    p6 <- p6 + geom_hline(yintercept = prec, color = "red", linetype = "dashed") + geom_vline(xintercept = rec, color = "red", linetype = "dashed")
    p6 <- p6 + geom_label(mapping = aes(label = paste0("-log10(pval) = ", round(-log10(cutoff), digits = 2))), hjust = "inward") + geom_point(color = "red")
    p6 <- p6 + geom_line(data = data.frame(precision = infsitesperf@y.values[[1]], recall = infsitesperf@x.values[[1]]), mapping = aes(x = recall, y = precision))
    p6 <- p6 + theme_minimal() + coord_equal(xlim = c(0,1), ylim = c(0,1)) + labs(y = "Precision", x = "Recall")
    # p6 <- p6 + labs(title = paste0("sample ", sampleid, " - ", "\ntotal testable: ", nrow(vafhitsdf), " - parallel violations: ", nrow(finalhits), " / total phaseable: ", nrow(vafhitsroc), " - violation confirmed: ", sum(vafhitsroc$is_confirmed)))
    
    ggsave(filename = file.path(sampledir, paste0(sampleid, "_PrecRec.png")), plot = p6, width = 10, height = 10, units = "cm")
  }
  # print(optimperf)
  return(optimperf)
}




# get_performance_metrics <- function(df, plotting = F, sampleid = sampleid, sampledir = OUTDIR) {
#   
#   if (sum(df$is_confirmed) == 0) {
#     outv <- setNames(object = numeric(length = 7L), nm = c("fdrcutoff", "fdrprec", "fdrrec", "Fone", "Fcutoff", "Fprec", "Frec"))
#     # print(outv)
#     return(outv)
#   }
#   
#   infsitespred <- prediction(predictions = -log10(df$pval), labels = df$is_confirmed)
#   infsitesperf <- performance(prediction.obj = infsitespred, measure = "prec", x.measure = "rec")
#   
#   fdridx <- tail(which(infsitesperf@y.values[[1]] >= .9), n = 1)
#   if (length(fdridx) == 0) {
#     fdridx <- which.max(infsitesperf@y.values[[1]])
#   }
#   fdrcutoff <- 10^-infsitesperf@alpha.values[[1]][fdridx]
#   fdrprec <- infsitesperf@y.values[[1]][fdridx]
#   fdrrec <- infsitesperf@x.values[[1]][fdridx]
#   
#   infsitesperfF <- performance(prediction.obj = infsitespred, measure = "f")
#   Fidx <- which.max(infsitesperfF@y.values[[1]])
#   Fone <- infsitesperfF@y.values[[1]][Fidx]
#   Fcutoff <- 10^-infsitesperfF@x.values[[1]][Fidx]
#   Fprec <- infsitesperf@y.values[[1]][Fidx]
#   Frec <- infsitesperf@x.values[[1]][Fidx]
#   
#   optimperf <- c(fdrcutoff = fdrcutoff, fdrprec = fdrprec, fdrrec = fdrrec, Fone = Fone, Fcutoff = Fcutoff, Fprec = Fprec, Frec = Frec)
#   
#   if (plotting) {
#     p6 <- ggplot(data = as.data.frame(t(optimperf)), mapping = aes(y = Fprec, x = Frec))
#     p6 <- p6 + geom_hline(yintercept = .9, color = "blue", linetype = "dashed") + geom_vline(xintercept = fdrrec, color = "blue", linetype = "dashed")
#     p6 <- p6 + geom_hline(yintercept = Fprec, color = "red", linetype = "dashed") + geom_vline(xintercept = Frec, color = "red", linetype = "dashed")
#     p6 <- p6 + geom_label(mapping = aes(label = paste0("-log10(pval) = ", round(-log10(Fcutoff), digits = 2))), nudge_x = .05, hjust = "left") + geom_point(color = "red")
#     p6 <- p6 + geom_label(mapping = aes(y = fdrprec, x = fdrrec, label = paste0("-log10(pval) = ", round(-log10(fdrcutoff), digits = 2))), nudge_x = .05, hjust = "left") + geom_point(mapping = aes(y = fdrprec, x = fdrrec), color = "blue")
#     p6 <- p6 + geom_line(data = data.frame(precision = infsitesperf@y.values[[1]], recall = infsitesperf@x.values[[1]]), mapping = aes(x = recall, y = precision))
#     p6 <- p6 + theme_minimal() + coord_equal(xlim = c(0,1), ylim = c(0,1)) + labs(y = "Precision", x = "Recall")
#     # p6 <- p6 + labs(title = paste0("sample ", sampleid, " - ", "\ntotal testable: ", nrow(vafhitsdf), " - parallel violations: ", nrow(finalhits), " / total phaseable: ", nrow(vafhitsroc), " - violation confirmed: ", sum(vafhitsroc$is_confirmed)))
#     
#     ggsave(filename = file.path(sampledir, paste0(sampleid, "_PrecRec.png")), plot = p6, width = 10, height = 10, units = "cm")
#   }
#   # print(optimperf)
#   return(optimperf)
# }
# 
