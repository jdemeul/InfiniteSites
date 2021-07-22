### new simulator for infinite sites violations in samples
### uses samples from same tumour type with same mut profile to resample half of the mut and determine violation counts

library(ggplot2)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(VariantAnnotation)

### start by reading in mutational profiles and determining which samples can be resampled from which
# OUTDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20200603_SampleBasedSim/"
OUTDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210617_SampleBasedSim_NODRIVERS/"
isaviolreffile <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isavioltypes_reffile.txt"
# dir.create(path = OUTDIR, showWarnings = F)

NSIMS <- 1000
NCORES <- 12





###### additional functions

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


### core function to be repeated
# sample mutload / 2 SNVs from similar samples and keep only the ones which are infsites violations
sample_single_iteration <- function(mutload, snvstosample, overlapsobj, violrefidx, bases, unused = NA) {
  
  infsitesidx <- intersect(x = sample(x = 1:length(snvstosample), size = mutload/2, replace = F), y = subjectHits(overlapsobj))
  insitesviolations <- snvstosample[infsitesidx[runif(n = length(infsitesidx)) > .5]] # every SNV has a prob of .5 of actually being picked from the original sample
  
  idxv <- setNames(object = 1:length(bases$trinucleotides_mutations), nm = bases$trinucleotides_mutations)
  mut1 <- idxv[insitesviolations$mut1]
  mut2 <- idxv[insitesviolations$mut2]
  alhitsdf <- data.frame(mut1 = bases$trinucleotides_mutations[pmin(mut1, mut2)], mut2 = bases$trinucleotides_mutations[pmax(mut1, mut2)])
  isaviolcounts <- data.frame(table(mut1 = factor(alhitsdf$mut1, levels = bases$trinucleotides_mutations), 
                                    mut2 = factor(alhitsdf$mut2, levels = bases$trinucleotides_mutations)))[violrefidx, "Freq"]

  return(isaviolcounts)
}






simulate_infsites_sample_based <- function(sampleid, allsnvmnv, useablesamples, nsims = NSIMS, outbase = OUTDIR, isaviolref, ncores = NCORES, bases) {
  
  # sampleid <- "9fc5b5c7-3973-42b4-8710-454de0cb5b50"
  sampleoutdir <- file.path(outbase, sampleid)
  dir.create(sampleoutdir, showWarnings = F)
  
  
  # get a tumour and its snvs and preprocess
  tumsnv <- allsnvmnv[[sampleid]]
  tumsnv$ALT <- unlist(tumsnv$ALT, recursive = T)
  strand(tumsnv) <- ifelse(tumsnv$REF %in% c("C", "T"), "+", "-")
  tumsnv$trinuc <- getSeq(resize(width = 3, x = tumsnv, fix = "center"), x = BSgenome.Hsapiens.1000genomes.hs37d5)
  tumsnv$mut1 <- paste0(tumsnv$trinuc, ">", ifelse(strand(tumsnv) == "+", tumsnv$ALT, reverseComplement(tumsnv$ALT)))
  
  # get all sampleable tumours and compile their snvs
  snvstosample <- unlist(allsnvmnv[useablesamples[[sampleid]]], use.names = F)
  snvstosample$ALT <- unlist(x = snvstosample$ALT, recursive = T)
  # snvstosample$REF <- as.character(x = snvstosample$REF)
  # snvstosample$ALT1 <- ""
  snvstosample$mut1 <- ""
  snvstosample$mut2 <- ""
  
  # compile table with all sampleable snvs with additional colum stating type of violation and its muttypes 
  tumor_sampleset_overlap <- findOverlaps(query = tumsnv, subject = snvstosample, type = "equal")
  
  # if there are no overlaps, no point in continuing
  if (length(tumor_sampleset_overlap) > 0 & length(snvstosample) >= length(tumsnv)) {
    
    snvstosample$mut1[subjectHits(tumor_sampleset_overlap)] <- tumsnv$mut1[queryHits(tumor_sampleset_overlap)]
    snvstosample$mut2[subjectHits(tumor_sampleset_overlap)] <- paste0(substr(x = snvstosample$mut1[subjectHits(tumor_sampleset_overlap)], start = 1, stop = 4), 
                                                                      ifelse(snvstosample$REF[subjectHits(tumor_sampleset_overlap)] %in% c("C", "T"), 
                                                                             snvstosample$ALT[subjectHits(tumor_sampleset_overlap)], 
                                                                             reverseComplement(snvstosample$ALT[subjectHits(tumor_sampleset_overlap)])))
    
    
    isaviolations <- mclapply(X = 1:nsims, FUN = sample_single_iteration, mutload = length(tumsnv), snvstosample = snvstosample, overlapsobj = tumor_sampleset_overlap, violrefidx = isaviolref$idx, bases = bases, mc.preschedule = T, mc.cores = ncores)
    
    finalcounts <- cbind(isaviolref[, -1], t(apply(X = do.call(cbind, isaviolations), MARGIN = 1, FUN = quantile, probs = c(.025, .5, .975), simplify = T)), Reduce('+', isaviolations))
    colnames(finalcounts) <- c("mut1", "mut2", "type", "freqal_low", "freqal_med", "freqal_hi", "counts_tot")
    
    
    # and per type
    isaviolations_pertype <- lapply(X = isaviolations, FUN = function(x, types) c(by(data = x, INDICES = types, FUN = sum)), types = isaviolref$type)
    isaviolations_pertype <- data.frame(t(apply(X = do.call(rbind, isaviolations_pertype), MARGIN = 2, FUN = quantile, probs = c(.025, .5, .975), simplify = T)))
    isaviolations_pertype$type <- row.names(isaviolations_pertype)
    colnames(isaviolations_pertype) <- c("freq_low", "freq_med", "freq_hi", "type")
    
  } else {
    finalcounts <- data.frame(isaviolref[, -1], freqal_low = 0, freqal_med = 0, freqal_hi = 0, counts_tot = 0)
    isaviolations_pertype <- data.frame(type = c("back", "forward", "parallel", "third_allele"), freq_low = 0, freq_med = 0, freq_hi = 0)
  }
  
  
  #writing output
  outfile1 <- file.path(sampleoutdir, paste0(sampleid, "_infsites_totals_samplebasedsim.txt"))
  outfile2 <- file.path(sampleoutdir, paste0(sampleid, "_infsites_permut_samplebasedsim.txt"))
  
  write.table(x = isaviolations_pertype[, c( "type", "freq_low", "freq_med", "freq_hi")], file = outfile1, quote = F, sep = "\t", row.names = F)  
  write.table(x = finalcounts, file = outfile2, quote = F, sep = "\t", row.names = F)
  
  return(NULL)
}





###### Setting up data and performing basic checks of feasibility
isaviolref <- read.delim(file = isaviolreffile, as.is = T)

RELEASETABLEFILE <- "/camp/project/proj-emedlab-vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
reltable <- read.delim(file = RELEASETABLEFILE, as.is = T)
combinedannots <- read.delim("/camp/project/proj-emedlab-vanloo/ICGC_annotations/summary_table_combined_annotations_v2.txt", as.is = T)
preferred_whitelisted_samples <- intersect(x =  combinedannots[combinedannots$is_preferred, "samplename"],
                                           y = unlist(strsplit(x = reltable[which(reltable$wgs_exclusion_white_gray == "Whitelist"), "tumor_wgs_aliquot_id"], split = ",")))



#### UPDATE 20210616 remove crappy PRAD-CA samples at this point
bad_PRAD_CA_samples <- c("08e1d976-6c39-428e-a4c2-f655b675683e", "108b67d4-5d66-46da-8675-6acae643b76f", "120f01d1-8884-4aca-a1cb-36b207b2aa3a", "126ee433-d345-4cac-882a-c91831a24690", "48c33a30-557b-4ecf-8066-5b4b068b5e3a", "4c5228b5-bf31-4abd-a47c-d088e16dba13", 
                         "4d11d7da-1204-437e-87b1-e8337a67c9a8", "61a48c69-4f7d-4dc6-aff7-88a6c33137df", "6d936ef9-b5df-44d3-831f-528bf8ddc131", "7ae9b843-488f-459c-8c0d-c81dcae57f99", "7dd2dc62-0eb4-4d45-86f1-e9e9377181ca", "887616c5-06a7-4e83-948c-3546202349fb", 
                         "a08ec059-7592-4698-bb45-25a9c3680c23", "ab8a55ed-ff47-4cad-ad91-52b9dc25aca7", "b2ec0fd0-fbcf-4abc-ad80-4ae444e30b55", "b33978c6-a855-4f9d-a0b0-d79453b9de41", "b33b7c8f-0b0d-4009-88a7-48e9d9cae6cb", "c08f65a0-bf4c-462e-9d07-ad56b3adcac8", 
                         "dcc938da-3e45-4c2f-ae0f-47817be04518", "e41bc2ec-3e0b-4c37-806b-3f6f25c8c4db", "f0f2030e-17fd-4dd9-9104-899e59d72ed8", "f5378545-17d4-4a64-a57e-f6c91ef4cb3a", "f640d377-98e9-41d3-8761-61eb33072c65")
preferred_whitelisted_samples <- setdiff(x = preferred_whitelisted_samples, y = bad_PRAD_CA_samples)
####



SIGINSAMPLESFILE <- "/camp/project/proj-emedlab-vanloo/ICGC_signatures/20180322_release/PCAWG_sigProfiler_SBS_signatures_in_samples_waliqID.csv"
sigs <- read.csv(file = SIGINSAMPLESFILE, as.is = T)
# drop all greylisted and excluded samples
sigs <- sigs[which(sigs$Sample.Name %in% preferred_whitelisted_samples), ]

# compute cosine sims
matprod <- as.matrix(sigs[, -c(1:3)]) %*% as.matrix(t(sigs[, -c(1:3)]))
vnorms <- sqrt(rowSums(sigs[, -c(1:3)]^2))
prodvnorms <- vnorms %*% t(vnorms)
cossim <- matprod/prodvnorms
rownames(cossim) <- sigs$Sample.Name
colnames(cossim) <- sigs$Sample.Name

# plot(1:2780, cossim["5c3def3a-b515-41f6-8157-681b963534e7",])
rm(matprod, vnorms, prodvnorms)

# for every sample, get vector of samples which can mutations can be sampled from 
minsimilarity <- 0.9
# similarsamples <- apply(X = cossim, MARGIN = 1, FUN = function(x, minsim) names(which(x >= minsim)), minsim = minsimilarity)

# plotdf <- data.frame(sigs[,1:2], sim = cossim["5c3def3a-b515-41f6-8157-681b963534e7",], row.names = NULL)
# p1 <- ggplot(data = plotdf, mapping = aes(x = 1:nrow(plotdf), y = sim, colour = Cancer.Types)) + geom_point()
# p1

## QC, check distribution of scores per disease type: e.g. mel vs lung

# for every sample, get other samples with similar mut spectrum of same cancer type
useablesamples <- lapply(X = sigs$Sample.Name, FUN = function(sid, sidctype, sims, minsim) {
  tumourtypesamples <- sidctype[sidctype$Cancer.Types == sidctype[sidctype$Sample.Name == sid, "Cancer.Types"], "Sample.Name"]
  simsigsamples <- names(which(sims[sid, ] >= minsim))
  outsamples <- setdiff(x = intersect(tumourtypesamples, simsigsamples), y = sid)
  return(outsamples)
}, sidctype = sigs[,1:2], sims = cossim, minsim = minsimilarity)

names(useablesamples) <- sigs$Sample.Name

mutloadpersample <- setNames(object = rowSums(x = sigs[, -c(1:3)]), nm = sigs$Sample.Name)

#generate overview DF for checking what will be possible
checkdf <- data.frame(cancertype = sigs$Cancer.Types, sampleid = sigs$Sample.Name, nsimsamples = lengths(useablesamples), mutload = mutloadpersample, samplespace = sapply(X = useablesamples, FUN = function(x, loads) sum(loads[x]), loads = mutloadpersample))
checkdf$nfold <- checkdf$samplespace / (checkdf$mutload/2)



### now that we have the tumours which we can simulate using this approach and which tumours to sample from, we can start loading the actual mutations

bases <- generate_bases_types_trinuc()



#### Addition 20210616 >> remove drivers
alldrivers <- read.delim(file = "/camp/project/proj-emedlab-vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
alldriversgr <- GRanges(seqnames = alldrivers$chr, ranges = IRanges(start = alldrivers$pos, width = nchar(alldrivers$ref)), REF = alldrivers$ref, ALT = alldrivers$alt, sampleid = alldrivers$sample, seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
alldriversgr <- alldriversgr[alldriversgr$REF %in% c("A", "C", "G", "T") & alldriversgr$ALT %in% c("A", "C", "G", "T")]
names(alldriversgr) <- paste0(alldriversgr$sampleid, "_", seqnames(alldriversgr), ":", start(alldriversgr), "_", alldriversgr$REF, "/", alldriversgr$ALT)

filterdrivers <- T
#####



# vcffiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/", pattern = ".somatic.snv_mnv.vcf.gz$", full.names = T, recursive = T)
vcffiles <- paste0("/camp/project/proj-emedlab-vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sigs$Sample.Name, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
all(file.exists(vcffiles))
# scanVcfHeader(vcffiles[1])
allsnvmnv <- GRangesList(mclapply(X = vcffiles, FUN = function(x) rowRanges(readVcf(file = x, param = ScanVcfParam(fixed = c("ALT"), info = NA, geno = NA))), mc.cores = NCORES, mc.preschedule = T))
names(allsnvmnv) <- sigs$Sample.Name

if (filterdrivers) {
  #### Addition 20210616 >> remove drivers
  allsnvmnv_nodriv <- GRangesList(mclapply(X = names(allsnvmnv), FUN = function(x, snvs, driv) {
    print(paste0("Running ", x))
    outgr <- subsetByOverlaps(x = snvs[[x]], ranges = driv[driv$sampleid %in% x], invert = T, type = "equal")
    return(outgr)
  }, snvs = allsnvmnv, driv = alldriversgr, mc.cores = NCORES, mc.preschedule = T))
  names(allsnvmnv_nodriv) <- names(allsnvmnv)
  allsnvmnv <- allsnvmnv_nodriv
  rm(allsnvmnv_nodriv, alldrivers, alldriversgr)
}



# some checks - found one snv-mnv file is missing, likely blacklisted at late stage, remove from useable samples
# notinsigfile <- setdiff(x = names(allsnvmnv), y = sigs$Sample.Name)
# nosnvfile <- setdiff(y = names(allsnvmnv), x = sigs$Sample.Name)
# useablesamples <- lapply(X = useablesamples, FUN = function(x, nofiles) x[!x %in% nofiles], nofiles = nosnvfile)
# 
# nosnvfile
#### this is where the fun stuff starts
## savving env image
save.image(file = file.path(OUTDIR, "allsnvs_NODRIVERS"))
# load(file = file.path(OUTDIR, "allsnvs_NODRIVERS.RDS"))


# debug(simulate_infsites_sample_based)
for (sampleid in checkdf[checkdf$nsimsamples > 0 & checkdf$nfold >= 1, "sampleid"]) {
  print(paste0("Running sample ", sampleid))
  simulate_infsites_sample_based(sampleid = sampleid, allsnvmnv = allsnvmnv, useablesamples = useablesamples, nsims = NSIMS, outbase = OUTDIR, isaviolref = isaviolref, ncores = NCORES, bases = bases)
}



