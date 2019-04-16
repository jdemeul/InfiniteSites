### functions for identifying gene conversion based on the allele frequencies

# get tumour allelecounts 
get_tum_allecounts_gr <- function(allelecountsdir, sampleid, bsgenome, reference_alleles, tempdir = TEMPDIR) {
  allelecountsfile_tum <- file.path(allelecountsdir, paste0(sampleid, "_allelecounts.tar.gz"))
  if (file.exists(allelecountsfile_tum)) {
    allelecounts_tum <- load_allelecounts(allelecountsfile = allelecountsfile_tum, scratchdir = tempdir)
  } else {
    return(NULL)
  }
  allelecounts_tum_gr <- GRanges(seqnames = allelecounts_tum$chr, IRanges(start = allelecounts_tum$pos, end = allelecounts_tum$pos), seqinfo = seqinfo(bsgenome))
  mcols(allelecounts_tum_gr) <- allelecounts_tum[, -c(1,2)]
  mcols(allelecounts_tum_gr)[, c("ref", "alt")] <- mcols(reference_alleles)

  
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


# get phased BAFs (for QC purposes, plotting, validation)
get_phased_BAF <- function(bafdir, sampleid, bsgenome, allelecounts) {
  phasedbaf_file <- file.path(bafdir, paste0(sampleid, ".BAFsegmented.txt"))
  if (file.exists(phasedbaf_file)) {
    phasedbaf <- read_tsv(file = phasedbaf_file, col_types = "cinnn")
  } else {
    return(NULL)
  }
  phasedbaf_gr <- GRanges(seqnames = phasedbaf$Chromosome, IRanges(start = phasedbaf$Position, end = phasedbaf$Position), seqinfo = seqinfo(bsgenome))
  mcols(phasedbaf_gr) <- phasedbaf[, -c(1:2)]
  
  locimatches <- findOverlaps(query = phasedbaf_gr, subject = allelecounts)
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
  gccorrlogr_file <- file.path(logrdir, paste0(sampleid, "_allelecounts"), paste0(sampleid, "_mutantLogR_gcCorrected.tab"))
  if (file.exists(gccorrlogr_file)) {
    gccorrlogr <- read_tsv(file = gccorrlogr_file, col_types = "cin")
  } else {
    return(NULL)
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
  
  write_tsv(x = as.data.frame(segments_gr), path = file.path(sampledir, paste0(sampleid, "_baf_logr_summarised_segments.txt")))
  
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
  som_vcf <- som_vcf[which(!is.na(info(som_vcf)$t_alt_count) & !is.na(info(som_vcf)$t_ref_count))]
  som <- rowRanges(som_vcf)
  mcols(som) <- cbind(mcols(som), info(som_vcf))
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



test_clean_sites <- function(sampledir, sampleid, segments_gr, phasedbaf_gr, snvs_vcf, bsgenome, logr, minsegmentSNPcount = 100, ncores = 12, presched = T, pseudocountrange = c(50, 1000), subsample_optim = T, recompute_pseudocoverage = T, immune_loci) {

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
  
  # filter SNPS where power is lacking: not very effective in PCAWG setting since virtually all sites are powered
  # if (filterpval < 1) {
  #   pfilter <- pmin(test_betabin_model(size = mcols(phasedbaf_gr)$refCount + mcols(phasedbaf_gr)$altCount,
  #                                      q = rep(0, length(phasedbaf_gr)),
  #                                      shape1 = mcols(phasedbaf_gr)$BAFseg*1000,
  #                                      shape2 = (1 - mcols(phasedbaf_gr)$BAFseg)*1000,
  #                                      idxs = clean_idx),
  #                   test_betabin_model(size = mcols(phasedbaf_gr)$refCount + mcols(phasedbaf_gr)$altCount,
  #                                      q = mcols(phasedbaf_gr)$refCount + mcols(phasedbaf_gr)$altCount,
  #                                      shape1 = mcols(phasedbaf_gr)$BAFseg*1000,
  #                                      shape2 = (1 - mcols(phasedbaf_gr)$BAFseg)*1000,
  #                                      idxs = clean_idx))
  #   
  #   clean_idx <- clean_idx[which(pfilter <= 1e-5)]
  # }
  # 
  # pval <- test_betabin_model(size = mcols(phasedbaf_gr)$refCount + mcols(phasedbaf_gr)$altCount,
  #                               q = mcols(phasedbaf_gr)$MajCount, 
  #                               shape1 = mcols(phasedbaf_gr)$BAFseg*1000, 
  #                               shape2 = (1 - mcols(phasedbaf_gr)$BAFseg)*1000,
  #                               idxs = clean_idx)
  # 
  # pphased <- data.frame(pval = pval, padj = p.adjust(pval, method = "hochberg"))
  
  # pphased$
  
  ####### flip only when BAF segment stricly greater than .5, otherwise inflating pvals
  # pvalre <- test_betabin_model(size = mcols(phasedbaf_gr)$refCount + mcols(phasedbaf_gr)$altCount,
  #                              q = pmax(mcols(phasedbaf_gr)$refCount, mcols(phasedbaf_gr)$altCount), 
  #                              shape1 = mcols(phasedbaf_gr)$BAFseg*1000, 
  #                              shape2 = (1 - mcols(phasedbaf_gr)$BAFseg)*1000,
  #                              idxs = clean_idx)
  
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
  pvalphased <- mcmapply(size = bbparams$size,
                     q = bbparams$majcount,
                     shape1 = bbparams$shape1*pseudocount_calibr,
                     shape2 = bbparams$shape2*pseudocount_calibr,
                     FUN = betabinom.test.ab,
                     MoreArgs = list(alternative = "two.sided"),
                     mc.preschedule = presched,
                     mc.cores = ncores)

  ggd.qqplot(sampleid = sampleid, sampledir = sampledir, suffix = "_phasedSNPs_QQplot", pvector = pvalphased)
  ## QC of the fit above, not used further
  # 
  # pvalre <- mcmapply(size = bbparams$size,
  #                    q = bbparams$maxcount,
  #                    shape1 = bbparams$shape1*pseudocount_calibr,
  #                    shape2 = bbparams$shape2*pseudocount_calibr,
  #                    FUN = betabinom.test.ab,
  #                    MoreArgs = list(alternative = "two.sided"),
  #                    mc.preschedule = presched,
  #                    mc.cores = ncores)
  # 
  # padjre <- p.adjust(pvalre, method = "fdr")
  # 
  # ggd.qqplot(sampleid = sampleid, sampledir = sampledir, suffix = "_dephasedSNPs_QQplot", pvector = pvalre)
  
  ######## somatic variants
  # tag immune regions
  mcols(snvs_vcf)$immune_locus <- snvs_vcf %within% immune_loci
  
  # subset SNVs to those on (long) segments (and PAR when appriopriate)
  snvseghits <- findOverlaps(query = snvs_vcf, subject = segments_gr)
  snvs_vcf <- snvs_vcf[queryHits(snvseghits)]
  
  pvalsnv <- mcmapply(size = mcols(snvs_vcf)$t_alt_count + mcols(snvs_vcf)$t_ref_count,
                     q = mcols(snvs_vcf)$t_alt_count,
                     shape1 = mcols(segments_gr)[subjectHits(snvseghits), "mu_vaf"]*pseudocount_calibr,
                     shape2 = (1-mcols(segments_gr)[subjectHits(snvseghits), "mu_vaf"])*pseudocount_calibr,
                     FUN = betabinom.test.ab,
                     MoreArgs = list(alternative = "greater"),
                     mc.preschedule = presched,
                     mc.cores = ncores)
  
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
  mcols(snvs_vcf) <- cbind(mcols(snvs_vcf), mcols(segments_gr)[subjectHits(snvseghits), c("total_cn", "major_cn", "minor_cn") ])
  
  # annotate hits with upstream and downstream BAF/LogR of SNPS
  # BAF
  preidxs <- precede(x = snvs_vcf, subject = phasedbaf_gr)
  postidxs <- follow(x = snvs_vcf, subject = phasedbaf_gr)
  mcols(snvs_vcf)$bafpos_pre <- start(phasedbaf_gr)[postidxs]
  mcols(snvs_vcf)$bafpos_post <- start(phasedbaf_gr)[preidxs]
  mcols(snvs_vcf)$bafpval_pre <- pvalphased[postidxs]
  mcols(snvs_vcf)$bafpval_post <- pvalphased[preidxs]
  
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

  outdf <- data.frame(chr = seqnames(snvs_vcf), start = start(snvs_vcf), end = end(snvs_vcf), ref = as.character(snvs_vcf$REF), alt = as.character(unlist(snvs_vcf$ALT)))
  outdf <- cbind(outdf, as.data.frame(mcols(snvs_vcf)[, c("VAF", "t_alt_count", "t_ref_count", "snv_near_indel", "Variant_Classification", "immune_locus", "pval",
                                                          "total_cn", "major_cn", "minor_cn", "bafpos_pre", "bafpos_post", "bafpval_pre", "bafpval_post", 
                                                          "logrpos_pre", "logrpos_post", "logrpval_pre", "logrpval_post")]))
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
                             p_empirical = numeric()), path = file.path(sampledir, paste0(sampleid, "_baflogr_gc_results.txt")))
    print(paste0("Sample ", sampleid, ": no conversion events detected"))
    return(hitvariants)
  }
  
  if (plotting)
    dir.create(file.path(sampledir, "figures"), showWarnings = F)
  
  hitvariants_grlist_annotated <- lapply(X = split(hitvariants, f = 1:length(hitvariants)), FUN = annotate_gc_hit, pseudocount_calibr = pseudocount_calibr, sampledir = sampledir, sampleid = sampleid, baf = baf, logr = logr, seg_gr = seg_gr, snvs = snvs, plotting = plotting, conversionlength = conversionlength, testwindow = testwindow, fdr = fdr)
  hitvariants_gr_annotated <- unlist(GRangesList(hitvariants_grlist_annotated))
  
  write_tsv(x = as.data.frame(hitvariants_gr_annotated), path = file.path(sampledir, paste0(sampleid, "_baflogr_gc_results.txt")))
  
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

