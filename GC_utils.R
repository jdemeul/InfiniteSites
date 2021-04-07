## wrapper for pbetabinomial ~ binom.test
betabinom.test <- function (q, size, prob, rho, alternative = c("two.sided", "less", "greater")) 
{
  # DNAME <- deparse(substitute(q))
  # qr <- round(q)
  # if (any(is.na(q) | (q < 0)) || max(abs(q - qr)) > 1e-07)
  #   stop("'q' must be nonnegative and integer")
  # q <- qr
  # if (length(q) == 2L) {
  #   size <- sum(q)
  #   q <- q[1L]
  # }
  # else if (length(q) == 1L) {
  #   sizer <- round(size)
  #   if ((length(size) > 1L) || is.na(size) || (size < 1) || abs(size - sizer) > 1e-07 || (q > sizer))
  #     stop("'size' must be a positive integer >= 'q'")
  #   DNAME <- paste(DNAME, "and", deparse(substitute(size)))
  #   size <- sizer
  # }
  # else stop("incorrect length of 'q'")
  # if (!missing(prob) && (length(prob) > 1L || is.na(prob) || prob < 0 ||
  #                     prob > 1))
  #   stop("'prob' must be a single number between 0 and 1")
  # alternative <- match.arg(alternative)
  PVAL <- switch(alternative, less = pbetabinom(q = q, size = size, prob = prob, rho = rho), greater = 1 - pbetabinom(q = q - 1, size = size, prob = prob, rho = rho),
                 two.sided = { if (prob == 0) (q == 0) else if (prob == 1) (q == size) else {
                   relErr <- 1 + 1e-07
                   d <- dbetabinom(x = q, size = size, prob = prob, rho = rho)
                   m <- size * prob
                   if (q == m) 1 else if (q < m) {
                     i <- seq.int(from = ceiling(m), to = size)
                     y <- sum(dbetabinom(x = i, size = size, prob = prob, rho = rho) <= d * relErr)
                     p1 <- 1 - pbetabinom(q = size - y, size = size, prob = prob, rho = rho)
                     if (p1 < 0)
                       pbetabinom(q = q, size = size, prob = prob, rho = rho)
                     else
                       p1 + pbetabinom(q = q, size = size, prob = prob, rho = rho) 
                   } else {
                     i <- seq.int(from = 0, to = floor(m))
                     y <- sum(dbetabinom(x = i, size = size, prob = prob, rho = rho) <= d * relErr)
                     p1 <- 1 - pbetabinom(q = q - 1, size = size, prob = prob, rho = rho)
                     if (p1 < 0)
                       pbetabinom(q = y - 1, size = size, prob = prob, rho = rho)
                     else 
                       p1 + pbetabinom(q = y - 1, size = size, prob = prob, rho = rho)
                     
                   }
                 }
                 })
}

## wrapper for pbetabinomial ~ binom.test
betabinom.test.ab <- function (q, size, shape1, shape2, alternative = c("two.sided", "less", "greater")) 
{
  DNAME <- deparse(substitute(q))
  qr <- round(q)
  sizer <- round(size)
  if (is.na(q) || (q < 0) || abs(q - qr) > 1e-07) {
    warning("'q' must be nonnegative and integer")
    return(NA)
  } else if ((length(size) > 1L) || is.na(size) || (size < 1) || abs(size - sizer) > 1e-07 || (q > sizer)) {
    warning("'size' must be a positive integer >= 'q'")
    return(NA)
  } else if (is.na(shape1) || is.na(shape2)){
    warning("'shapeX' must be nonnegative integers")
    return(NA)
  }
  q <- qr
  size <- sizer
  # if (length(q) == 2L) {
  #   size <- sum(q)
  #   q <- q[1L]
  # }
  # else if (length(q) == 1L) {
  #   if ((length(size) > 1L) || is.na(size) || (size < 1) || abs(size - sizer) > 1e-07 || (q > sizer))
  #     stop("'size' must be a positive integer >= 'q'")
  #   DNAME <- paste(DNAME, "and", deparse(substitute(size)))
  # }
  # else stop("incorrect length of 'q'")
  # if (!missing(prob) && (length(prob) > 1L || is.na(prob) || prob < 0 ||
  #                     prob > 1))
  #   stop("'prob' must be a single number between 0 and 1")
  
  alternative <- match.arg(alternative)
  PVAL <- switch(alternative, less = pbetabinom.ab(q = q, size = size, shape1 = shape1, shape2 = shape2), greater = 1 - pbetabinom.ab(q = q - 1, size = size, shape1 = shape1, shape2 = shape2),
                 two.sided = { if (shape1 == 0) (q == 0) else if (shape2 == 0) (q == size) else {
                   relErr <- 1 + 1e-07
                   d <- dbetabinom.ab(x = q, size = size, shape1 = shape1, shape2 = shape2)
                   m <- size * shape1 / (shape1 + shape2)
                   if (q == m) 1 else if (q < m) {
                     i <- seq.int(from = ceiling(m), to = size)
                     y <- sum(dbetabinom.ab(x = i, size = size, shape1 = shape1, shape2 = shape2) <= d * relErr)
                     p1 <- 1 - pbetabinom.ab(q = size - y, size = size, shape1 = shape1, shape2 = shape2)
                     if (p1 < 0)
                       pbetabinom.ab(q = q, size = size, shape1 = shape1, shape2 = shape2)
                     else
                       p1 + pbetabinom.ab(q = q, size = size, shape1 = shape1, shape2 = shape2) 
                   } else {
                     i <- seq.int(from = 0, to = floor(m))
                     y <- sum(dbetabinom.ab(x = i, size = size, shape1 = shape1, shape2 = shape2) <= d * relErr)
                     p1 <- 1 - pbetabinom.ab(q = q - 1, size = size, shape1 = shape1, shape2 = shape2)
                     if (p1 < 0)
                       pbetabinom.ab(q = y - 1, size = size, shape1 = shape1, shape2 = shape2)
                     else 
                       p1 + pbetabinom.ab(q = y - 1, size = size, shape1 = shape1, shape2 = shape2)
                     
                   }
                 }
                 })
}


test_betabin_model <- function(size, q, shape1, shape2, idxs, alternative = "two.sided", presched = T, ncores = NCORES) {
  pval <- mcmapply(size = size[idxs], q = q[idxs], shape1 = shape1[idxs], shape2 = shape2[idxs],
                   FUN = betabinom.test.ab, MoreArgs = list(alternative = alternative), mc.preschedule = presched, mc.cores = ncores)
  return(pval)
}
# 
# test_betabin_model_filter <- function(size, shape1, shape2, idxs, alternative = "two.sided", presched = T, ncores = NCORES) {
#   pval <- mcmapply(size = size[idxs], q = size[idxs], shape1 = shape1[idxs], shape2 = shape2[idxs],
#                    FUN = betabinom.test.ab, MoreArgs = list(alternative = alternative), mc.preschedule = presched, mc.cores = ncores)
#   return(pval)
# }


ggd.qqplot = function(sampleid = "", sampledir = ".", suffix = "_QQplot", pvector, main=NULL, tofile = F, ...) {
  o = -log10(sort(pvector,decreasing=F))
  o <- o[is.finite(o)]
  e = -log10( 1:length(o)/length(o) )
  png(filename = file.path(sampledir, paste0(sampleid, suffix, ".png")))
  plot(e,o,pch=".",cex=2, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
  dev.off()
}


ggd_qqplot2 = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  o <- o[is.finite(o)]
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=".",cex=2, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}


qqplot2 = function(pvector, pvector_conserv, main=NULL, ...) {
  ord <- order(pvector, decreasing = F)
  o = -log10(pvector[ord])
  o_cons <- -log10(pvector_conserv[ord])
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=".",cex=2, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  points(e, o_cons, col="blue", cex=1)
  lines(e,e,col="red")
}

read_pcawg_release_table <- function(release_table_file) {
  release_table <- read.delim(file = release_table_file, as.is = T)
  release_table <- release_table[release_table$wgs_exclusion_white_gray != "Excluded", ]
  splitaliquots <- strsplit(x = release_table$tumor_wgs_aliquot_id, split = ",")
  release_table_dedup <- release_table[rep(1:nrow(release_table), lengths(splitaliquots)), c("wgs_exclusion_white_gray", "dcc_project_code", "normal_wgs_aliquot_id", "tumor_wgs_aliquot_id", "normal_wgs_bwa_alignment_bam_file_name", "tumor_wgs_bwa_alignment_bam_file_name")]
  release_table_dedup$tumor_wgs_aliquot_id <- unlist(splitaliquots)
  release_table_dedup$tumor_wgs_bwa_alignment_bam_file_name <- unlist(strsplit(x = release_table$tumor_wgs_bwa_alignment_bam_file_name, split = ","))
  return(release_table_dedup)
}



# difference betwene two betas P(B2 > B1), see Chris Stucchio

h <- function(alpha_1, beta_1, alpha_2, beta_2) {
  j <- seq.int(0, round(alpha_2)-1)
  log_vals <- (lbeta(alpha_1 + j, beta_1 + beta_2) - log(beta_2 + j) -
                 lbeta(1 + j, beta_2) - lbeta(alpha_1, beta_1))
  1 - sum(exp(log_vals))
}



betaABfromMeanKappa = function( mean , kappa ) {
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}


betaABfromMeanSD = function( mean , sd ) {
  kappa = mean*(1-mean)/sd^2 - 1
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( c( a=a , b=b ) )
}


load_1000G_reference_alleles <- function(refallelesdir, chrominfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)) {
  refallelefiles <- paste0(refallelesdir, "1000genomesAlleles2012_chr", 1:23, ".txt")
  refalleles <- lapply(X = refallelefiles, function(x) read_tsv(file = x, col_types = "iii"))
  chr <- rep(c(1:22, "X"), sapply(X = refalleles, FUN = nrow))
  refalleles <- as.data.frame(do.call(what = rbind, args = refalleles))
  refalleles_gr <- GRanges(seqnames = chr, IRanges(start = refalleles$position, end = refalleles$position), seqinfo = chrominfo)
  mcols(refalleles_gr)$ref <- factor(refalleles$a0, levels = 1:4, labels = c("A", "C", "G", "T"))
  mcols(refalleles_gr)$alt <- factor(refalleles$a1, levels = 1:4, labels = c("A", "C", "G", "T"))
  return(refalleles_gr)
}


HDIofICDF = function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {
  # Arguments:
  #   ICDFname is R's name for the inverse cumulative density function
  #     of the distribution.
  #   credMass is the desired mass of the HDI region.
  #   tol is passed to R's optimize function.
  # Return value:
  #   Highest density iterval (HDI) limits in a vector.
  # Example of use: For determining HDI of a beta(30,12) distribution, type
  #   HDIofICDF( qbeta , shape1 = 30 , shape2 = 12 )
  #   Notice that the parameters of the ICDFname must be explicitly named;
  #   e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
  # Adapted and corrected from Greg Snow's TeachingDemos package.
  incredMass =  1.0 - credMass
  intervalWidth = function( lowTailPr , ICDFname , credMass , ... ) {
    ICDFname( credMass + lowTailPr , ... ) - ICDFname( lowTailPr , ... )
  }
  optInfo = optimize( intervalWidth , c( 0 , incredMass ) , ICDFname=ICDFname ,
                      credMass=credMass , tol=tol , ... )
  HDIlowTailPr = optInfo$minimum
  return( c( ICDFname( HDIlowTailPr , ... ) ,
             ICDFname( credMass + HDIlowTailPr , ... ) ) )
}


get_mean_var <- function(x) {
  return(c(mu = mean(x, na.rm = T), rho = sqrt(var(x, na.rm = T)) ))
}

get_mean_var_robust <- function(x) {
  return(c(mu = median(x, na.rm = T), rho = mad(x, na.rm = T) ))
}


