### new simulator for infinite sites violations in samples
### uses samples from same tumour type with same mut profile to resample half of the mut and determine violation counts

library(ggplot2)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(VariantAnnotation)

### start by reading in mutational profiles and determining which samples can be resampled from which
OUTDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210616_SampleBasedSim_pertumourtype_NODRIVERS/"
isaviolreffile <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isavioltypes_reffile.txt"
# dir.create(path = OUTDIR, showWarnings = F)

NSIMS <- 1000
NCORES <- 8




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
sample_single_iteration_ttype <- function(mutload, allsnvmnv, overlapsobj, bases, violrefidx, nsims, idxv) {
  # this is wehere a simulation starts
  sampledidx <- sample(x = 1:length(allsnvmnv), size = mutload, replace = F)
  
  infsiteshits <- overlapsobj[queryHits(overlapsobj) %in% sampledidx & subjectHits(overlapsobj) %in% sampledidx]
  # get rid of triple hits and above
  infsiteshits <- infsiteshits[!duplicated(queryHits(infsiteshits))]
  
  mut1 <- idxv[allsnvmnv$mut[queryHits(infsiteshits)]]
  mut2 <- idxv[allsnvmnv$mut[subjectHits(infsiteshits)]]
  alhitsdf <- data.frame(mut1 = bases$trinucleotides_mutations[pmin(mut1, mut2)], mut2 = bases$trinucleotides_mutations[pmax(mut1, mut2)])
  isaviolcounts <- data.frame(table(mut1 = factor(alhitsdf$mut1, levels = bases$trinucleotides_mutations), mut2 = factor(alhitsdf$mut2, levels = bases$trinucleotides_mutations)))[violrefidx, "Freq"]
  
  return(isaviolcounts)
}



simulate_infsites_ttype <- function(ttype, allsnvmnv, nsims = NSIMS, outbase = OUTDIR, isaviolref, ncores = NCORES, bases, mutloads) {
  
  sampleoutdir <- file.path(outbase, ttype)
  dir.create(sampleoutdir, showWarnings = F)
  
  # statics
  idxv <- setNames(object = 1:length(bases$trinucleotides_mutations), nm = bases$trinucleotides_mutations)
  
  # preprocess
  allsnvmnv$ALT <- unlist(allsnvmnv$ALT, recursive = T)
  strand(allsnvmnv) <- ifelse(allsnvmnv$REF %in% c("C", "T"), "+", "-")
  allsnvmnv$trinuc <- getSeq(resize(width = 3, x = allsnvmnv, fix = "center"), x = BSgenome.Hsapiens.1000genomes.hs37d5)
  allsnvmnv$mut <- paste0(allsnvmnv$trinuc, ">", ifelse(strand(allsnvmnv) == "+", allsnvmnv$ALT, reverseComplement(allsnvmnv$ALT)))
  
  # compile table with all sampleable snvs with additional colum stating type of violation and its muttypes 
  tumor_sampleset_overlap <- findOverlaps(query = allsnvmnv, type = "equal", drop.self = T, drop.redundant = T)
  
  # debug(sample_single_iteration_ttype)
  # allfinalcounts <- list()
  # allisaviolations_pertype <- list()
  for (mutload in mutloads) {
    if (mutload > length(allsnvmnv)/2) next
    
    isaviolations <- mclapply(X = 1:nsims, FUN = sample_single_iteration_ttype, mutload = mutload, allsnvmnv = allsnvmnv, overlapsobj = tumor_sampleset_overlap, violrefidx = isaviolref$idx, bases = bases, idxv = idxv, mc.preschedule = T, mc.cores = ncores)
    
    
    finalcounts <- cbind(isaviolref[, -1], t(apply(X = do.call(cbind, isaviolations), MARGIN = 1, FUN = quantile, probs = c(.025, .5, .975), simplify = T)), Reduce('+', isaviolations))
    colnames(finalcounts) <- c("mut1", "mut2", "type", "freqal_low", "freqal_med", "freqal_hi", "counts_tot")
    
    # and per type
    isaviolations_pertype <- lapply(X = isaviolations, FUN = function(x, types) c(by(data = x, INDICES = types, FUN = sum)), types = isaviolref$type)
    isaviolations_pertype <- data.frame(t(apply(X = do.call(rbind, isaviolations_pertype), MARGIN = 2, FUN = quantile, probs = c(.025, .5, .975), simplify = T)))
    isaviolations_pertype$type <- row.names(isaviolations_pertype)
    colnames(isaviolations_pertype) <- c("freq_low", "freq_med", "freq_hi", "type")
    
    # allfinalcounts[[paste0("load", mutload)]] <- finalcounts
    # allisaviolations_pertype[[paste0("load", mutload)]] <- isaviolations_pertype
    
    #writing output
    outfile1 <- file.path(sampleoutdir, paste0(ttype, "_load_", mutload, "_infsites_totals_samplebasedsim.txt"))
    outfile2 <- file.path(sampleoutdir, paste0(ttype, "_load_", mutload, "_infsites_permut_samplebasedsim.txt"))
    
    write.table(x = isaviolations_pertype[, c( "type", "freq_low", "freq_med", "freq_hi")], file = outfile1, quote = F, sep = "\t", row.names = F)  
    write.table(x = finalcounts, file = outfile2, quote = F, sep = "\t", row.names = F)
  }
  
  
  
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

## QC, check distribution of scores per disease type: e.g. mel vs lung

# for every sample, get other samples with similar mut spectrum of same cancer type
useablesamples <- lapply(X = sigs$Sample.Name, FUN = function(sid, sidctype, sims, minsim) {
  tumourtypesamples <- sidctype[sidctype$Cancer.Types == sidctype[sidctype$Sample.Name == sid, "Cancer.Types"], "Sample.Name"]
  simsigsamples <- names(which(sims[sid, ] >= minsim))
  outsamples <- intersect(tumourtypesamples, simsigsamples)
  return(outsamples)
}, sidctype = sigs[,1:2], sims = cossim, minsim = minsimilarity)

names(useablesamples) <- sigs$Sample.Name

mutloadpersample <- setNames(object = rowSums(x = sigs[, -c(1:3)]), nm = sigs$Sample.Name)

#generate overview DF for checking what will be possible
checkdf <- data.frame(cancertype = sigs$Cancer.Types, sampleid = sigs$Sample.Name, nsimsamples = lengths(useablesamples), mutload = mutloadpersample, samplespace = sapply(X = useablesamples, FUN = function(x, loads) sum(loads[x]), loads = mutloadpersample), stringsAsFactors = F)
checkdf$nfold <- checkdf$samplespace / (checkdf$mutload/2)

representativedf <- data.frame(do.call(rbind, by(data = checkdf, INDICES = checkdf$cancertype, FUN = function(x) x[which.max(x$nsimsamples), ])), stringsAsFactors = F) # may want to add cluster of the ultrahypermuts
representativedf <- representativedf[representativedf$samplespace >= 1e5 & representativedf$nsimsamples >= 10,]
write.table(x = representativedf, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210616_SampleBasedSim_pertumourtype_representativesamples", sep = "\t", quote = F, row.names = F)


### now that we have the tumours which we can simulate using this approach and which tumours to sample from, we can start loading the actual mutations

bases <- generate_bases_types_trinuc()



#### Addition 20210616 >> remove drivers
alldrivers <- read.delim(file = "/camp/project/proj-emedlab-vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
alldriversgr <- GRanges(seqnames = alldrivers$chr, ranges = IRanges(start = alldrivers$pos, width = nchar(alldrivers$ref)), REF = alldrivers$ref, ALT = alldrivers$alt, sampleid = alldrivers$sample, seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
alldriversgr <- alldriversgr[alldriversgr$REF %in% c("A", "C", "G", "T") & alldriversgr$ALT %in% c("A", "C", "G", "T")]
names(alldriversgr) <- paste0(alldriversgr$sampleid, "_", seqnames(alldriversgr), ":", start(alldriversgr), "_", alldriversgr$REF, "/", alldriversgr$ALT)


filterdrivers <- T



#### main function
# debug(simulate_infsites_sample_based)
for (ttype in representativedf$cancertype) {
  # ttype <- "Uterus-AdenoCA"
  
  sampleid <- representativedf[representativedf$cancertype == ttype, "sampleid"]
  # vcffiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/", pattern = ".somatic.snv_mnv.vcf.gz$", full.names = T, recursive = T)
  vcffiles <- paste0("/camp/project/proj-emedlab-vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", useablesamples[[sampleid]], ".consensus.20160830.somatic.snv_mnv.vcf.gz")
  allsnvmnv_init <- mclapply(X = vcffiles, FUN = function(x) rowRanges(readVcf(file = x, param = ScanVcfParam(fixed = c("ALT"), info = NA, geno = NA))), mc.cores = 6, mc.preschedule = T)
  allsnvmnv_init <- unlist(GRangesList(mapply(allsnvmnv_init, useablesamples[[sampleid]], FUN = function(x, y) {x$sampleid <- y; return(x)})))
  # allsnvmnv_init <- unlist(GRangesList(mclapply(X = vcffiles, FUN = function(x) rowRanges(readVcf(file = x, param = ScanVcfParam(fixed = c("ALT"), info = NA, geno = NA))), mc.cores = 6, mc.preschedule = T)))
  
  if (filterdrivers) {
    # relevantdrivers <- alldriversgr[alldriversgr$sampleid %in% sampleid]
    # findOverlap(query = allsnvmnv_init, subject = )
    idx <- paste0(allsnvmnv_init$sampleid, "_", seqnames(allsnvmnv_init), ":", start(allsnvmnv_init), "_", allsnvmnv_init$REF, "/", unlist(allsnvmnv_init$ALT))
    allsnvmnv_init <- allsnvmnv_init[!idx %in% names(alldriversgr[alldriversgr$sampleid %in% useablesamples[[sampleid]]])]
  }
  
  
  print(paste0("Running ", ttype))
  # debug(simulate_infsites_ttype)
  simulate_infsites_ttype(ttype = ttype, allsnvmnv = allsnvmnv_init, mutloads = 2^c(10:20), nsims = NSIMS, outbase = OUTDIR, isaviolref = isaviolref, ncores = NCORES, bases = bases)
}
#######




####### downstream analyses and QC

### quick sanity checks
# load new simualtor results
# newsimres <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20200603_SampleBasedSim/", pattern = "_infsites_permut_samplebasedsim.txt", full.names = T, recursive = T)
newsimres <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210617_SampleBasedSim_NODRIVERS/", pattern = "_infsites_permut_samplebasedsim.txt", full.names = T, recursive = T)
newsimresls <- mclapply(X = newsimres, FUN = read.delim, as.is = T, mc.cores = 6, mc.preschedule = T)
names(newsimresls) <- gsub(pattern = "_infsites_permut_samplebasedsim.txt", replacement = "", x = basename(newsimres))

originalres <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_all_writefrac/", pattern = "_infsites_permut_effgenfrac1.txt", full.names = T, recursive = T)
# originalres <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_het_only_writefrac/", pattern = "_infsites_permut_effgenfrac0.1.txt", full.names = T, recursive = T)
originalresls <- mclapply(X = originalres, FUN = read.delim, as.is = T, mc.cores = 6, mc.preschedule = T)
names(originalresls) <- gsub(pattern = "_infsites_permut_effgenfrac1.txt", replacement = "", x = basename(originalres))

jointsamples <- intersect(names(newsimresls), names(originalresls))

outcosine <- setNames(object = rep(0, length(jointsamples)), nm = jointsamples)
outtotals <- setNames(object = rep(0, length(jointsamples)), nm = jointsamples)
outtotals_new <- setNames(object = rep(0, length(jointsamples)), nm = jointsamples)
bialidx <- which(originalresls[[1]]$type %in% c("parallel", "third_allele"))
for (sampleid in jointsamples) {
  # sampleid <- jointsamples[1]
  vorig <- as.numeric(originalresls[[sampleid]]$counts_tot[bialidx])
  vnew <- as.numeric(newsimresls[[sampleid]]$counts_tot[bialidx])
  outcosine[[sampleid]] <- sum(vnew*vorig) / ( sqrt(sum(vorig^2)) * sqrt(sum(vnew^2)) )
  outtotals[[sampleid]] <- sum(vorig)
  outtotals_new[[sampleid]] <- sum(vnew)
}

hist(outcosine[outtotals > 100], breaks = 40)
df1 <- data.frame(sampleid = names(outtotals), total = outtotals, total_new = outtotals_new, cosine = outcosine)
df1 <- merge(x = df1, y = sigs[, 1:2], by.x = "sampleid", by.y = "Sample.Name")


source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(df1$Cancer.Types))), scheme = "tumour.subtype")
names(cvect) <- unique(df1$Cancer.Types)
# cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Skin-Melanoma-Mucosal", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", "#000000", '#FF4500', '#FF4500')

bad_PRAD_CA_samples <- c("08e1d976-6c39-428e-a4c2-f655b675683e", "108b67d4-5d66-46da-8675-6acae643b76f", "120f01d1-8884-4aca-a1cb-36b207b2aa3a", "126ee433-d345-4cac-882a-c91831a24690", "48c33a30-557b-4ecf-8066-5b4b068b5e3a", "4c5228b5-bf31-4abd-a47c-d088e16dba13", 
                         "4d11d7da-1204-437e-87b1-e8337a67c9a8", "61a48c69-4f7d-4dc6-aff7-88a6c33137df", "6d936ef9-b5df-44d3-831f-528bf8ddc131", "7ae9b843-488f-459c-8c0d-c81dcae57f99", "7dd2dc62-0eb4-4d45-86f1-e9e9377181ca", "887616c5-06a7-4e83-948c-3546202349fb", 
                         "a08ec059-7592-4698-bb45-25a9c3680c23", "ab8a55ed-ff47-4cad-ad91-52b9dc25aca7", "b2ec0fd0-fbcf-4abc-ad80-4ae444e30b55", "b33978c6-a855-4f9d-a0b0-d79453b9de41", "b33b7c8f-0b0d-4009-88a7-48e9d9cae6cb", "c08f65a0-bf4c-462e-9d07-ad56b3adcac8", 
                         "dcc938da-3e45-4c2f-ae0f-47817be04518", "e41bc2ec-3e0b-4c37-806b-3f6f25c8c4db", "f0f2030e-17fd-4dd9-9104-899e59d72ed8", "f5378545-17d4-4a64-a57e-f6c91ef4cb3a", "f640d377-98e9-41d3-8761-61eb33072c65")
origdf1 <- df1
df1 <- df1[!df1$sampleid %in% bad_PRAD_CA_samples, ]
p1 <- ggplot(data = df1, mapping = aes(x = total/1000, y = total_new/1000, size = cosine))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_point(shape = 21, mapping = aes(fill = Cancer.Types), color = "black", alpha = .4) + 
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + annotation_logticks() + scale_fill_manual(values = cvect, name = "Cancer type") + 
  labs(x = "Average biallelic violations per 1000 simulations (Permutation model)", y = "Average biallelic violations per 1000 simulations (Resampling model)", size = "Cosine similarity violation spectra") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4))
p1 <- p1 + theme(panel.grid = element_blank())
p1
ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_noPRAD_CA_NODRIV.pdf", width = 10, height = 8)
write.table(x = df1, file = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA_NODRIV.txt", quote = F, row.names = F, sep = "\t")
# ggsave(plot = p1, filename = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_noPRAD_CA.pdf", width = 10, height = 8)
# write.table(x = df1, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA.txt", quote = F, row.names = F, sep = "\t")
# all(df1$sampleid %in% preferred_whitelisted_samples)


### load one time only
callablegenome <- import.wig("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/data/pancan_callable_genome.wig")
callablegenomesize <- sum(width(callablegenome[mcols(callablegenome)$score == 1]))



# df1 <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA.txt", as.is = T)
# additional plot with boxplots to show per tumour type effects of ratio
ctypenamestouse <- names(which(table(df1[df1$total >= 10 & df1$total_new >= 10, "Cancer.Types"]) >= 10))
p2 <- ggplot(data = df1[df1$total > 0 & df1$total_new > 0 & df1$Cancer.Types %in% ctypenamestouse, ], mapping = aes(x = Cancer.Types, y = log10(callablegenomesize/(total_new/total)))) + 
  geom_hline(yintercept = log10(callablegenomesize), linetype = "dashed") + 
  geom_point(mapping = aes(colour = Cancer.Types), position = position_jitter(width = 0.25), shape = 19, show.legend = F, size = 0.5) + geom_boxplot(mapping = aes(fill = Cancer.Types), outlier.shape = NA, alpha = 0.3, show.legend = F) + labs(y = "Effective genome size (bp)") 
p2 <- p2 + scale_color_manual(values = cvect, name = "Cancer type") + scale_fill_manual(values = cvect, name = "Cancer type") + theme_minimal() + theme(axis.text = element_text(angle = 90), axis.title.y.right = element_text(angle = 90), axis.title.x = element_blank(), panel.grid = element_blank()) 
p2 <- p2 + annotation_logticks(sides = "r", scaled = T, base = 10) + scale_x_discrete(limits = sort(unique(ctypenamestouse))) + scale_y_continuous(breaks = c(7,8,9,10), position = "right") 
p2

# p2 <- ggplot(data = df1[df1$total > 0 & df1$total_new > 0 & df1$Cancer.Types %in% ctypenamestouse, ], mapping = aes(x = Cancer.Types, y = log2(total_new/total))) + 
#   geom_hline(yintercept = 0) +
#   geom_point(mapping = aes(colour = Cancer.Types), position = position_jitter(width = 0.25), shape = 19, show.legend = F, size = 0.5) + geom_boxplot(mapping = aes(fill = Cancer.Types), outlier.shape = NA, alpha = 0.3, show.legend = F) + labs(y = "log2(# resampling-based violations / # permutation-based violations)") 
# p2 <- p2 + scale_color_manual(values = cvect, name = "Cancer type") + scale_fill_manual(values = cvect, name = "Cancer type") + theme_minimal() + theme(axis.text.x = element_text(angle = 90), axis.title.y = element_blank(), panel.grid = element_blank()) 
# p2 <- p2 + scale_x_discrete(limits = rev(sort(unique(ctypenamestouse)))) + scale_y_continuous(breaks = c(-2,0,2,4,6,8)) + coord_flip()
# p2
ggsave(plot = p2, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_noPRAD_CA_boxplot_redo_NODRIV_revorder.pdf", width = 6, height = 4, useDingbats=FALSE)


df1check <- df1[df1$total > 0 & df1$total_new > 0 & df1$Cancer.Types %in% ctypenamestouse, ]
sizevect <- c(by(data = df1check, INDICES = df1check$Cancer.Types, FUN = function(x) summary(x$total_new/x$total)))
callablegenomesize/sizevect[["Lymph-BNHL"]]

compdf2 <- cbind(originalresls[["303abbe5-4155-4a0d-bc3b-f8995261ca52"]], newsimresls[["303abbe5-4155-4a0d-bc3b-f8995261ca52"]])




### note, a small set of ColoRect/Uterine samples seem to persistently have LOW cosine sims of their allelic infinite sites hits between simulators: these have in common a SBS44 exposure 
useablesamples[["303abbe5-4155-4a0d-bc3b-f8995261ca52"]]
temp <- c("303abbe5-4155-4a0d-bc3b-f8995261ca52", useablesamples[["303abbe5-4155-4a0d-bc3b-f8995261ca52"]])
sigs[sigs$Sample.Name %in% temp, ]

vcffiles <- paste0("/camp/project/proj-emedlab-vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", temp, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
file.exists(vcffiles)
sbs44vcfs <- lapply(X = paste0("/camp/project/proj-emedlab-vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", temp, ".consensus.20160830.somatic.snv_mnv.vcf.gz"), FUN = function(x) rowRanges(readVcf(file = x, param = ScanVcfParam(fixed = c("ALT"), info = NA, geno = NA))))
### majority of these are greylisted, due to lot of germline SNP calls (likely contamination)
RELEASETABLEFILE <- "/camp/project/proj-emedlab-vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
reltable <- read.delim(file = RELEASETABLEFILE, as.is = T)

df1$greylist <- reltable$wgs_exclusion_white_gray[match(x = df1$sampleid, table = reltable$tumor_wgs_aliquot_id)]

df1$greylist[is.na(df1$greylist)] <- "excluded"

p1 <- ggplot(data = df1[df1$total > 1000, ], mapping = aes(x = greylist, y = cosine)) + geom_violin() + geom_jitter(mapping = aes(color = Cancer.Types, alpha = log10(total+1)), width = 0.1)
p1

p1 <- ggplot(data = df1, mapping = aes(x = log10(total+1), y = cosine)) + geom_point(mapping = aes(color = Cancer.Types, shape = greylist))
p1

# lesson learned: do NOT use even graylisted samples for this





#### quick plot of curve

#### quick plot of curve

# melresfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210616_SampleBasedSim_pertumourtype/", pattern = "_infsites_totals_samplebasedsim.txt", full.names = T, recursive = T)
# melresfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210616_SampleBasedSim_pertumourtype/", pattern = "_infsites_permut_samplebasedsim.txt", full.names = T, recursive = T)
melresfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210616_SampleBasedSim_pertumourtype_NODRIVERS/", pattern = "_infsites_permut_samplebasedsim.txt", full.names = T, recursive = T)
# temp <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210616_SampleBasedSim_pertumourtype_NODRIVERS/Biliary-AdenoCA/Biliary-AdenoCA_load_65536_infsites_permut_samplebasedsim.txt", as.is = T)

melres <- lapply(X = melresfiles, FUN = function(x) {
  dfin <- read.delim(x, as.is = T)
  # df <- as.data.frame(t(colSums(df[3:4,2:4])))
  fileinfo <- strsplit(x = basename(x), split = "_")[[1]]
  # df[, c("tumour_histology", "simulated_burden")] <- fileinfo[c(1,3)]
  df <- data.frame(freq = sum(dfin$counts_tot)/1000, tumour_histology = fileinfo[1], simulated_burden = fileinfo[3])
  return(df)
})

melresdf <- do.call(rbind, melres)
melresdf$simulated_burden <- as.integer(melresdf$simulated_burden)


library(rtracklayer)
library(ggplot2)
callablegenome <- import.wig("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/data/pancan_callable_genome.wig")
callablegenomesize <- sum(width(callablegenome[mcols(callablegenome)$score == 1]))

source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(melresdf$tumour_histology))), scheme = "tumour.subtype")
names(cvect) <- unique(melresdf$tumour_histology)

fitcoefficients <- c(by(data = melresdf[melresdf$freq > 0.01, ], INDICES = melresdf[melresdf$freq > 0.01, "tumour_histology"], FUN = function(x) mean(x$freq / (x$simulated_burden^2 / callablegenomesize))))

# fitcoefficients <- do.call(rbind, by(data = melresdf, INDICES = melresdf$tumour_histology, FUN = function(x) x[which.max(x$simulated_burden), ]))
# fitcoefficients$fit <- fitcoefficients$freq/(fitcoefficients$simulated_burden^2/callablegenomesize)

options(scipen = 999)
p1 <- ggplot(data = rbind(melresdf, data.frame(freq = NA, tumour_histology = "ColoRect-AdenoCA", simulated_burden = 2*1048576)), mapping = aes(x = simulated_burden)) +
  geom_point(mapping = aes(y = freq, colour = tumour_histology), alpha = .6, show.legend = F, size = 2) +
  stat_function(fun = function(x) x^2/callablegenomesize, color = "black", linetype = "dashed", alpha = .8, )
p1 <- p1 + scale_x_log10(labels = 10^(4:6), breaks = 10^(4:6)) + scale_y_log10(labels = 10^(0:4), breaks = 10^(0:4))
p1 <- p1 + theme_minimal() + theme(panel.grid.minor = element_blank()) + annotation_logticks() + scale_colour_manual(values = cvect) + coord_cartesian(xlim = c(10^3, 1.5*10^6), ylim = c(10^-2, 2*10^4))
p1 <- p1 + labs(x = "Total SNV burden", y = "Biallelic mutation frequency")
# p1 <- p1 + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3))


p1 <- p1 + stat_function(fun = function(x) fitcoefficients[1]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[1]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[2]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[2]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[3]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[3]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[4]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[4]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[5]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[5]])

p1 <- p1 + stat_function(fun = function(x) fitcoefficients[6]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[6]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[7]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[7]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[8]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[8]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[9]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[9]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[10]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[10]])

p1 <- p1 + stat_function(fun = function(x) fitcoefficients[11]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[11]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[12]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[12]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[13]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[13]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[14]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[14]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[15]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[15]])

p1 <- p1 + stat_function(fun = function(x) fitcoefficients[16]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[16]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[17]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[17]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[18]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[18]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[19]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[19]])
p1 <- p1 + stat_function(fun = function(x) fitcoefficients[20]*x^2/callablegenomesize, color = cvect[names(fitcoefficients)[20]])

# p1 <- p1 + geom_point(data = df1, mapping = aes(x = mutburden, y = ifelse(tot_phaseable>1e4, med, nparallel) + ndivergent, colour = histology_abbreviation))

p1

ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Lineplot_coefficients_n2overN_NODRIVERS.pdf", width = 3, height = 6)



RELEASETABLEFILE <- "/camp/project/proj-emedlab-vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
reltable <- read.delim(file = RELEASETABLEFILE, as.is = T)

repsamples <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210616_SampleBasedSim_pertumourtype_representativesamples", as.is = T)
repsamples$icgc_donor_id <- reltable[match(x = repsamples$sampleid, table = reltable$tumor_wgs_aliquot_id), "icgc_donor_id"]

fitcoefficientsdf <- data.frame(coefficient = fitcoefficients, tumour_histology = names(fitcoefficients))
fitcoefficientsdf$icgc_donor_id <- repsamples[match(x = fitcoefficientsdf$tumour_histology, table = repsamples$cancertype), "icgc_donor_id"] 
fitcoefficientsdf$xlab <- paste0(fitcoefficientsdf$tumour_histology, " (", fitcoefficientsdf$icgc_donor_id, ")")
fitcoefficientsdf$xlab <- factor(x = fitcoefficientsdf$xlab, levels = fitcoefficientsdf$xlab[order(fitcoefficients, decreasing = T)])

p1 <- ggplot(data = fitcoefficientsdf, mapping = aes(x = xlab, y = log2(coefficient), fill = tumour_histology)) + geom_col(show.legend = F) + 
  scale_fill_manual(values = cvect) + theme_minimal() + theme(axis.text.x = element_text(angle = 90), panel.grid.major.x = element_blank())
p1 <- p1 + geom_text(mapping = aes(label = round(coefficient, digits = 1)), nudge_y = 0.1, angle = 90)
p1 <- p1 + geom_hline(yintercept = 0, linetype = "dashed") + labs(y = "Coefficient of n2/N", x = "") + 
  # scale_y_log10(labels = c(1,2,4,8,16,32,64), breaks = c(1,2,4,8,16,32,64), minor_breaks = c(1:64)) +
  annotation_logticks(base = 2, sides = "l") + 
  scale_y_continuous(labels = c(1,2,4,8,16,32,64), breaks = log2(c(1,2,4,8,16,32,64)), minor_breaks = NULL)
p1

ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/barplot_coefficients_n2overN_NODRIVERS.pdf", width = 6, height = 3)


repsamples$coefficient <- fitcoefficients[repsamples$cancertype]
write.table(x = repsamples, file = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/NeighbResampling_TumourType_plus_coefficients_n2overN_NODRIVERS.txt", quote = F, row.names = F, sep = "\t")




