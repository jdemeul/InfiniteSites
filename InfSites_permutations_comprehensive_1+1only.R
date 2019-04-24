## novel permutation strategy

# start with actual sample
# sample mut load (rgamma + rpois with mut load as lambda)
# sample mut load per chrom (dirichlet + rmultinom)
# sample mut types on chrom (dirichlet + rmultinom)


library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
# library(ggplot2)
library(VariantAnnotation)
library(parallel)
library(rslurm)
library(VGAM)
# library(ggplot2)

# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")

SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/"
RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
OUTDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_1plus1_only/"
CNDIR <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"

NCORES <- 12
NSIMS <- 100


gammaShRaFromModeSD = function( mode , sd ) {
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
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


get_mutspectrum <- function(mutations, trinucs) {
  
  agidxs <- which(as.character(mutations$REF) %in% c("A", "G"))
  mcols(mutations)[agidxs, "trinuc"] <- reverseComplement(mcols(mutations)[agidxs, "trinuc"])
  mcols(mutations)[agidxs, "ALT"] <- reverseComplement(mcols(mutations)[agidxs, "ALT"])
  mcols(mutations)[agidxs, "REF"] <- reverseComplement(mcols(mutations)[agidxs, "REF"])
  
  mcols(mutations)$type <- factor(paste0(mcols(mutations)$trinuc, ">", mcols(mutations)$ALT),
                                  levels = trinucs$trinucleotides_mutations)
  
  muttypefreqperchr <- table(seqnames(mutations), mcols(mutations)$type)
  
  return(muttypefreqperchr)
}


# classify/count isaviolations
classify_violations <- function(sampledsnvs) {
  
  # back and forward mutation
  bfhits <- findOverlaps(sampledsnvs, type = "equal", drop.self = T, drop.redundant = T)
  if (length(bfhits) > 0) {
    mut1 <- factor(x = paste0(sampledsnvs[queryHits(bfhits)]$trinucobs, ">", sampledsnvs[queryHits(bfhits)]$alt), levels = trinucs$trinucleotides_mutations)
    mut2 <- factor(x = paste0(sampledsnvs[subjectHits(bfhits)]$trinucobs, ">", sampledsnvs[subjectHits(bfhits)]$alt), levels = trinucs$trinucleotides_mutations)
    bfhitsdf <- data.frame(table(mut1, mut2))
    colnames(bfhitsdf) <- c("mut1", "mut2", "freqbf")
  } else {
    bfhitsdf <- expand.grid(mut1 = factor(trinucs$trinucleotides_mutations), mut2 = factor(trinucs$trinucleotides_mutations), freqbf = 0)
  }
  
  # hits to another allele
  sampleranges_nostrand <- subsetByOverlaps(x = sampledsnvs, ranges = sampledsnvs[queryHits(bfhits)], type = "equal", invert = T)
  strand(sampleranges_nostrand) <- "*"
  alhits <- findOverlaps(sampleranges_nostrand, type = "equal", drop.self = T, drop.redundant = T)
  if (length(alhits) > 0) {
    idxv <- setNames(object = 1:length(trinucs$trinucleotides_mutations), nm = trinucs$trinucleotides_mutations)
    mut1 <- idxv[paste0(sampleranges_nostrand[queryHits(alhits)]$trinucobs, ">", sampleranges_nostrand[queryHits(alhits)]$alt)]
    mut2 <- idxv[paste0(sampleranges_nostrand[subjectHits(alhits)]$trinucobs, ">", sampleranges_nostrand[subjectHits(alhits)]$alt)]
    alhitsdf <- data.frame(mut1 = trinucs$trinucleotides_mutations[pmin(mut1, mut2)], mut2 = trinucs$trinucleotides_mutations[pmax(mut1, mut2)])
    alhitsdf <- data.frame(table(mut1 = factor(alhitsdf$mut1, levels = trinucs$trinucleotides_mutations), mut2 = factor(alhitsdf$mut2, levels = trinucs$trinucleotides_mutations)))
    colnames(alhitsdf) <- c("mut1", "mut2", "freqal")
  } else {
    alhitsdf <- expand.grid(mut1 = factor(trinucs$trinucleotides_mutations), mut2 = factor(trinucs$trinucleotides_mutations), freqal = 0)
  }
  
  isaviolations <- cbind(bfhitsdf, freqal = alhitsdf$freqal)
  
  return(isaviolations)
}




simulate_snvs_chr <- function(chridx, nmutsim_chr, consecutivesegments, trinucs, mutspectrum_chr_norm, sizefactor = 25) {
  
  # print(paste0("chr", chridx))
  
  # get probabilities for types and acceptance rate for mutations
  probspertrinuc <- c(by(data = unlist(mutspectrum_chr_norm[chridx,]), INDICES = trinucs$trinucleotides, FUN = sum))
  acceptpertrinuc <- probspertrinuc/max(probspertrinuc)
  
  probsmuttype <- split(x = setNames(object = unlist(mutspectrum_chr_norm[chridx,]), nm = substr(x = colnames(mutspectrum_chr_norm), 5,5)), f = trinucs$trinucleotides)
  probsmuttype <- lapply(probsmuttype, FUN = function(x) x/sum(x))
  
  
  inithitidxs <- integer()
  while (length(inithitidxs) < 1.01*nmutsim_chr[[chridx]]) {
    
    # sample mappable positions
    sampleranges <- GRanges(seqnames = chridx, ranges = IRanges(start = sample(x = 1:max(end(consecutivesegments[[chridx]])), size = sizefactor*nmutsim_chr[[chridx]], replace = T), width = 3), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
    # shift to original segments
    sampleranges <- shift(x = sampleranges, shift = mcols(consecutivesegments[[chridx]])[nearest(x = sampleranges, subject = consecutivesegments[[chridx]]), "offset"] - 1)
    
    # sampleranges <- GRanges(seqnames = chridx, ranges = IRanges(start = sample(x = 1:seqlengths(BSgenome.Hsapiens.1000genomes.hs37d5)[[chridx]], size = sizefactor*nmutsim_chr[[chridx]], replace = T), width = 3), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
    # sampleranges <- subsetByOverlaps(x = sampleranges, ranges = callablegenome)
    
    # if (length(sampleranges) == 0) {
    #   print("Tricky region, sampling again")
    #   next
    # }
    
    # keep track of "time"
    sampleranges$timestep <- 1:length(sampleranges)
    
    # annotate with trinuc context and copy (encoded as strand here)
    sampleranges$trinuc <- getSeq(sampleranges, x = BSgenome.Hsapiens.1000genomes.hs37d5)
    agidxs <- which(as.character(subseq(sampleranges$trinuc, 2, 2)) %in% c("A", "G"))
    mcols(sampleranges)[agidxs, "trinuc"] <- reverseComplement(mcols(sampleranges)[agidxs, "trinuc"])
    strand(sampleranges) <- sample(x = c("+", "-"), size = length(sampleranges), replace = T)
    
    # omit trinucs containing N's
    naidxs <- grep(pattern = "N", x = sampleranges$trinuc)
    if (length(naidxs) > 0) {
      sampleranges <- sampleranges[-naidxs]
    }
    
    # uniform sampling to check for acceptance of mutation
    sampleranges$initunif <- runif(n = length(sampleranges), min = 0, max = 1)
    
    # extend with new iteration and shift timestep
    if (exists(x = "sampleranges_previousit")) {
      sampleranges$timestep <- sampleranges$timestep + length(sampleranges_previousit)
      sampleranges <- c(sampleranges_previousit, sampleranges)
    }
    
    inithitidxs <- which(sampleranges$initunif < acceptpertrinuc[as.character(sampleranges$trinuc)])
    
    sampleranges_previousit <- sampleranges
  }
  
  # identify all double hits (i.e. same position & same copy)
  doublehits <- findOverlaps(query = sampleranges, type = "equal", drop.self = T, drop.redundant = T)
  
  # keep only those which are initially hit OR a potential double hit (meaning, first time around the position was hit, reevaluate second time)
  inithits <- sampleranges[inithitidxs]
  potdoublehits <- sampleranges[unique(subjectHits(doublehits)[which(queryHits(doublehits) %in% inithitidxs)])]
  # remove any second hit from the initial hits for re-evaluation as second hit to new trinuc context
  inithits <- inithits[which(!inithits$timestep %in% potdoublehits$timestep)]
  
  # cleanup
  rm(sampleranges)
  
  # annotate mutations accepted during first pass
  inithits$alt <- sapply(X = as.character(inithits$trinuc), FUN = function(x, mutlist) sample(x = names(mutlist[[x]]), size = 1, prob = mutlist[[x]]), mutlist = probsmuttype)
  inithits$trinucalt <- replaceAt(x = inithits$trinuc, at = IRanges(start = 2, width = 1), value = CharacterList(split(inithits$alt, f = 1:length(inithits))))
  agidxs <- which(inithits$alt %in% c("A", "G"))
  mcols(inithits)[agidxs, "trinucalt"] <- reverseComplement(mcols(inithits)[agidxs, "trinucalt"])
  inithits$trinucobs <- inithits$trinuc
  mcols(inithits) <- mcols(inithits)[, c("timestep", "trinuc", "trinucobs", "alt", "trinucalt")]
  
  # match second hits to initial ones and set reference trinuc to alt for potential second hits
  matchedhits <- findOverlaps(query = potdoublehits, subject = inithits, type = "equal")
  potdoublehits$trinucobs <- mcols(inithits)[subjectHits(matchedhits), "trinucalt"]
  potdoublehits <- potdoublehits[which(potdoublehits$initunif < acceptpertrinuc[as.character(potdoublehits$trinucobs)])]
  
  if (length(potdoublehits) > 0) {
    # annotate actual second hits
    potdoublehits$alt <- sapply(X = as.character(potdoublehits$trinucobs), FUN = function(x, mutlist) sample(x = names(mutlist[[x]]), size = 1, prob = mutlist[[x]]), mutlist = probsmuttype)
    potdoublehits$trinucalt <- replaceAt(x = potdoublehits$trinucobs, at = IRanges(start = 2, width = 1), value = CharacterList(split(potdoublehits$alt, f = 1:length(potdoublehits))))
    agidxs <- which(potdoublehits$alt %in% c("A", "G"))
    mcols(potdoublehits)[agidxs, "trinucalt"] <- reverseComplement(mcols(potdoublehits)[agidxs, "trinucalt"])
    mcols(potdoublehits) <- mcols(potdoublehits)[, c("timestep", "trinuc", "trinucobs", "alt", "trinucalt")]
    
    # merging
    sampledsnvs <- c(inithits, potdoublehits)
  } else {
    # if no second hits set sampled snvs to the initial hits
    sampledsnvs <- inithits
  }
  
  sampledsnvs <- sampledsnvs[order(sampledsnvs$timestep, decreasing = F)]

    # get final set of mutations
    sampledsnvs <- sampledsnvs[1:nmutsim_chr[[chridx]]]
    
    # final classification
    isaviolations <- classify_violations(sampledsnvs = sampledsnvs)

  return(isaviolations)
}



simulate_snvs_sample_single_iteration <- function(mutspectrum, trinucs, callablegenome, ncores, callableseqtrinucsbychrom) {
  
  #### sampling repeats start here
  # resample total number of muts again
  gammashaperate <- gammaShRaFromModeSD(mode = sum(mutspectrum), sd = sum(mutspectrum)*.05)
  nmutsim <- rpois(n = 1, lambda = rgamma(n = 1, shape = gammashaperate$shape, rate = gammashaperate$rate))
  
  # resample number of muts per chrom
  nmutsim_chr <- c(rmultinom(n = 1, size = nmutsim, prob = MCMCpack::rdirichlet(n = 1, alpha = rowSums(mutspectrum)+1)))
  names(nmutsim_chr) <- rownames(mutspectrum)
  
  # omit chromosomes with ≤ 1 total mutations
  nmutsim_chr <- nmutsim_chr[nmutsim_chr > 1]
  
  # resample mut types per chrom
  muttypeprior <- (colSums(mutspectrum)+1)*ncol(mutspectrum)/sum(colSums(mutspectrum))
  mutspectrum_chr <- t(sapply(X = rownames(mutspectrum), FUN = function(chr, prior, obs) MCMCpack::rdirichlet(n = 1, alpha = mutspectrum[chr,]+prior),
                              obs = mutspectrum, prior = muttypeprior))
  colnames(mutspectrum_chr) <- colnames(mutspectrum)
  
  mutspectrum_chr_norm <- mutspectrum_chr/callableseqtrinucsbychrom[rownames(mutspectrum_chr), trinucs$trinucleotides, drop = F]
  mutspectrum_chr_norm <- mutspectrum_chr_norm/rowSums(mutspectrum_chr_norm)
  colnames(mutspectrum_chr_norm) <- trinucs$trinucleotides_mutations
  
  # # pre-generate samplevectors
  consecutivesegments <- lapply(X = names(nmutsim_chr), FUN = function(chridx, cn) {
    cn_chr <- cn[which(seqnames(cn) == chridx)]
    consecutivesegments <- GRanges(seqnames = chridx, ranges = IRanges(start = cumsum(c(1, width(cn_chr)[-length(cn_chr)])), width = width(cn_chr)), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
    mcols(consecutivesegments)$offset <- start(cn_chr) - start(consecutivesegments)
  return(consecutivesegments)
  }, cn = callablegenome)
  
  names(consecutivesegments) <- names(nmutsim_chr)

  
  #### per chromosome
  isaviolations_single_iteration <- mclapply(X = names(nmutsim_chr), FUN = simulate_snvs_chr, nmutsim_chr = nmutsim_chr, trinucs = trinucs, mutspectrum_chr_norm = mutspectrum_chr_norm, consecutivesegments = consecutivesegments, mc.preschedule = F, mc.cores = ncores)
  
  # combine and return
  isaviolations_single_iteration <- Reduce('+', lapply(isaviolations_single_iteration, FUN = function(x) x[, c("freqbf", "freqal")]))
  
  return(isaviolations_single_iteration)
}




generate_isa_breakdown_file <- function(sampleid, outdir, snvindeldir, cndir, nsims, ncores, callablegenome, trinucs) {

  print(paste0("Running sample ", sampleid))

  # load snvs
  snv_mnvfile <- list.files(path = snvindeldir, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz$"), full.names = T, recursive = T)
  cnfile <- paste0(cndir, sampleid, ".consensus.20170119.somatic.cna.annotated.txt")
  
  if (length(snv_mnvfile) == 0)
    return(NULL)
  
  snvs <- rowRanges(readVcf(file = snv_mnvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
  totalsnvs <- length(snvs)
  # always drop Y, keep X, also in males
  snvs <- snvs[lengths(mcols(snvs)$ALT) == 1 & seqnames(snvs) != "Y"]
  # MODS for 1+1
  cn <- GRanges(read.delim(file = cnfile, as.is = T), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
  cn <- reduce(cn[which(cn$major_cn == 1 & cn$minor_cn == 1)])
  snvs <- subsetByOverlaps(x = snvs, ranges = cn, type = "within")
  
  if (length(snvs) < 2)
    return(NULL)
  
  mcols(snvs)$ALT <- unlist(mcols(snvs)$ALT)
  
  mcols(snvs)$trinuc <- getSeq(x = BSgenome.Hsapiens.1000genomes.hs37d5, resize(snvs, width = 3, fix = "center"))
  
  mutspectrum <- get_mutspectrum(mutations = snvs, trinucs = trinucs)
  
  #
  callablegenomecn <- intersect(x = callablegenome, y = cn)
  callableseq <- getSeq(x = BSgenome.Hsapiens.1000genomes.hs37d5, callablegenomecn)
  callableseqtrinucsbychrom <- by(data = trinucleotideFrequency(callableseq), INDICES = seqnames(callablegenomecn), FUN = colSums)
  # callableseqtrinucs <- Reduce(f = '+', x = callableseqtrinucsbychrom)
  callableseqtrinucsbychrom <- data.frame(do.call(rbind, lapply(callableseqtrinucsbychrom, function(x) x/sum(x))))
  
  callableseqtrinucsbychrom <- callableseqtrinucsbychrom[, unique(trinucs$trinucleotides)] + 
    callableseqtrinucsbychrom[, as.character(reverseComplement(DNAStringSet(unique(trinucs$trinucleotides))))]
  #
  
  statsdf <- c(total_genome = sum(as.numeric(width(callablegenomecn))), frac_call_genome = sum(as.numeric(width(callablegenomecn))) / sum(as.numeric(width(callablegenome))), used_snvs = length(snvs), frac_snvs = length(snvs)/totalsnvs)
  write.table(x = statsdf, file = file.path(outdir, paste0(sampleid, "_generalstats.txt")), quote = F, sep = "\t", row.names = T, col.names = F)
  
  isaviolations <- mclapply(X = 1:nsims, FUN = function(x) {print(paste0("simulation ", x)); simulate_snvs_sample_single_iteration(mutspectrum = mutspectrum, trinucs = trinucs, callablegenome = callablegenomecn, callableseqtrinucsbychrom = callableseqtrinucsbychrom, ncores = 1)}, mc.allow.recursive = F, mc.preschedule = T, mc.cores = ncores)
  # isaviolations <- lapply(X = 1:nsims, FUN = function(x) {print(paste0("simulation ", x)); simulate_snvs_sample_single_iteration(mutspectrum = mutspectrum, trinucs = trinucs, callablegenome = callablegenome, callableseqtrinucsbychrom = callableseqtrinucsbychrom, ncores = 1)})
  # browser()
  bfcounts <- t(apply(X = do.call(cbind, lapply(isaviolations, function(x) x$freqbf)), MARGIN = 1, FUN = quantile, probs = c(.025, .5, .975), simplify = T))
  alcounts <- t(apply(X = do.call(cbind, lapply(isaviolations, function(x) x$freqal)), MARGIN = 1, FUN = quantile, probs = c(.025, .5, .975), simplify = T))
  finalcounts <- cbind(expand.grid(factor(trinucs$trinucleotides_mutations), factor(trinucs$trinucleotides_mutations)), bfcounts, alcounts)
  colnames(finalcounts) <- c("mut1", "mut2", "freqbf_low", "freqbf_med", "freqbf_hi", "freqal_low", "freqal_med", "freqal_hi")
  
  bftotals <- quantile(sapply(isaviolations, function(x) sum(x$freqbf)), probs = c(.025, .5, .975))
  altotals <- quantile(sapply(isaviolations, function(x) sum(x$freqal)), probs = c(.025, .5, .975))
  totals <- data.frame(mut1 = NA, mut2 = NA, t(bftotals), t(altotals))
  colnames(totals) <- c("mut1", "mut2", "freqbf_low", "freqbf_med", "freqbf_hi", "freqal_low", "freqal_med", "freqal_hi")
  
  finalcounts <- rbind(finalcounts, totals)
  
  write.table(x = finalcounts, file = file.path(outdir, paste0(sampleid, "_infsites_backfwd_allelic.txt")), quote = F, sep = "\t", row.names = F)
  saveRDS(object = isaviolations, file = file.path(outdir, paste0(sampleid, "_infsites_violations.rds")))
  
  return(finalcounts)
}





# sample positions on chrom (sample excess!)
# accept or reject candidate based on probs of trinuc (sum of all types from that trinuc) (normalised frequencies)
# sample muttype according to probs (i.e. normalised frequencies)
# when pos within trinuc of previous one, flag up and treat separately
#     same pos, 


### load one time only
callablegenome <- import.wig("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/data/pancan_callable_genome.wig")
seqlevelsStyle(callablegenome) <- "Ensembl"
callablegenome <- GRanges(callablegenome, seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
callablegenome <- sort(callablegenome[mcols(callablegenome)$score == 1])

trinucs <- generate_bases_types_trinuc()

# rslurmdf <- read.delim(SUMTABLE_WHOLE, as.is = T)
rslurmdf <- read_pcawg_release_table(release_table_file = RELEASETABLEFILE)[, "tumor_wgs_aliquot_id" , drop = F]
colnames(rslurmdf) <- "sampleid"

# do not redo stuff for now
# doneids <- gsub(pattern = "_infsites_backfwd_allelic.txt", replacement = "", x = list.files(path = OUTDIR, pattern = "_infsites_backfwd_allelic.txt$"))
# rslurmdf <- rslurmdf[!rslurmdf$sampleid %in% doneids,, drop = F]

# rm(callableseq, callableseqtrinucs)


# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
# sampleid <- "fc876f5c-8339-bc9c-e040-11ac0d485160"
# sampleid <- "7cae6c0b-36fe-411b-bbba-093a4c846d84"
# sampleid <- "0332b017-17d5-4083-8fc4-9d6f8fdbbbde"
# sampleid <- "f82d213f-bc06-5b51-e040-11ac0c48687e"


# debug(generate_isa_breakdown_file)
# undebug(simulate_snvs_chr)
# debug(simulate_snvs_sample_single_iteration)

# set.seed(seed = 13454)
# estisaviolations <- generate_isa_breakdown_file(sampleid = sampleid, outdir = OUTDIR, snvindeldir = SNVMNVINDELDIR, nsims = NSIMS,
#                             ncores = NCORES, callablegenome = callablegenome, trinucs = trinucs, cndir = CNDIR)

#### slurm wrapping
generate_isa_breakdown_file_slurmwrap <- function(sampleid) {
  estisaviolations <- generate_isa_breakdown_file(sampleid = sampleid, outdir = OUTDIR, snvindeldir = SNVMNVINDELDIR, cndir = CNDIR, nsims = NSIMS,
                                                  ncores = NCORES, callablegenome = callablegenome, trinucs = trinucs)
  return(NULL)
}

# debug(get_mutspectrum)
# generate_isa_breakdown_file_slurmwrap(sampleid = sampleid)

# 
isajob <- slurm_apply(f = generate_isa_breakdown_file_slurmwrap, params = rslurmdf[,,drop=F], jobname = "isajob4", nodes = 470, cpus_per_node = 1, add_objects = ls(),
                          pkgs = rev(.packages()), libPaths = .libPaths(), submit = T, slurm_options = list(exclude="fat-worker00[1-2]"))






