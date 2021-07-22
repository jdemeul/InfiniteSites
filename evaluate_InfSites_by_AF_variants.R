#### evaluate parallel calls

library(ggplot2)

source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")
CLEANHISTOLOGYFILE <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

# allsumfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", pattern = "_InfSites_VAFpipeline_summarystats.txt", full.names = T, recursive = T)
# allsumfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", pattern = "_InfSites_VAFpipeline_summarystats.txt", full.names = T, recursive = T)
allsumfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210702_vafpipeline_out_alphapt1_hetonly/", pattern = "_InfSites_VAFpipeline_summarystats.txt", full.names = T, recursive = T)
allsums <- lapply(X = allsumfiles, FUN = read.table, row.names = 1, as.is = T)
allsums <- data.frame(t(do.call(cbind, allsums)))
allsums$sampleid <- gsub(pattern = "_InfSites_VAFpipeline_summarystats.txt", replacement = "", x = basename(allsumfiles))

# allpcfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
# allpcfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
allpcfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210702_vafpipeline_out_alphapt1_hetonly/", pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
allpc <- lapply(X = allpcfiles, FUN = read.table, col.names = c("pseudocount", "snv_slope"), as.is = T)
allpc <- data.frame(do.call(rbind, allpc))
allpc$sampleid <- gsub(pattern = "_pseudocount_calibrated.txt", replacement = "", x = basename(allpcfiles))

allsums <- merge(x = allsums, y = allpc, by = "sampleid")

sumtab <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)
allsums$histology_abbreviation <- sumtab[match(x = allsums$sampleid, table = sumtab$samplename), "histology_abbreviation"]
consrhopsi <- read.delim(file = "/camp/project/proj-emedlab-vanloo/ICGC_consensus_copynumber/consensus.20170217.purity.ploidy.txt.gz", as.is = T)
allsums <- cbind(allsums, consrhopsi[match(x = allsums$sampleid, table = consrhopsi$samplename), -1])
allsums$is_preferred <- allsums$sampleid %in% sumtab[sumtab$is_preferred, "samplename"]


allsimfiles <- list.files(path = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_all_writefrac/", pattern = "_infsites_permut_effgenfrac1.txt", full.names = T, recursive = T)
allsims <- lapply(X = allsimfiles, FUN = function(x) {
  y <- read.table(x, as.is = T, header = T)
  nparsim <- sum(y[y$type == "parallel", "counts_tot"])/1000
  ndivsim <- sum(y[y$type == "third_allele", "counts_tot"])/1000
  return(c(nparsim = nparsim, ndivsim = ndivsim))
  })
allsims <- data.frame(do.call(rbind, allsims))
allsims$sampleid <- gsub(pattern = "_infsites_permut_effgenfrac1.txt", replacement = "", x = basename(allsimfiles))

allsums$simnumber <- allsims[match(x = allsums$sampleid, table = allsims$sampleid), "nparsim"]
allsums$simnumber_divergent <- allsims[match(x = allsums$sampleid, table = allsims$sampleid), "ndivsim"]
allsums$snvs_total <- rowSums(sumtab[match(x = allsums$sampleid, table = sumtab$samplename), c("num_clonal", "num_subclonal")])

# allsums$hitratio <- log10(allsums$nparallel/allsums$tot_testable)
# allsums$hitratio2 <- log10(allsums$med/allsums$tot_testable)
# allsums$hitratio3 <- log10(allsums$tot_testable^2/allsums$nparallel)
# allsums$hitratio4 <- log10(allsums$tot_testable^2/allsums$med)
# allsums$hitratio5 <- log10(allsums$tot_diploid^2/allsums$med_diploid)
# allsums$pratio <- pnorm(q = allsums$hitratio3, mean = mean(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), sd = sd(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), lower.tail = TRUE, log.p = FALSE)
# allsums$simratio <- log2((allsums$nparallel+0.001)/(allsums$simnumber+0.001))
# allsums$simratio2 <- log2((allsums$nparallel+0.001)/(allsums$simnumber*allsums$tot_testable/allsums$snvs_total+0.001))


# temp <- read.delim("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_het_only_writefrac/0009b464-b376-4fbc-8a56-da538269a02f/0009b464-b376-4fbc-8a56-da538269a02f_infsites_permut_effgenfrac0.1.txt", as.is = T)
# temp <- sum(temp[temp$type == "parallel", "counts_tot"])/1000
# p1 <- ggplot(data = allsums[which(allsums$prec > 0 & allsums$rec > 0), ]) + geom_point(mapping = aes(x = hitratio2, y = purity_conf_mad))
# p1

# p1 <- ggplot(data = allsums) + geom_point(mapping = aes(x = hitratio, y = purity_conf_mad, colour = hitratio > -2.352241))
# p1

# ratecutoff <- qnorm(p = .01, mean = mean(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), sd = sd(allsums[allsums$tot_phaseable >= 1000 & allsums$med > 0, "hitratio4"]), lower.tail = TRUE, log.p = FALSE)
# p1 <- ggplot(data = allsums) + geom_point(mapping = aes(y = log10(nparallel), x = log10(snvs_total), color = histology_abbreviation, shape = log2((nparallel+0.001)/(simnumber+0.001)) < 8), size = 5) +
#   annotation_logticks() + scale_color_manual(values = cvect)
# p1


### note
RELEASETABLEFILE <- "/camp/project/proj-emedlab-vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
reltable <- read.delim(file = RELEASETABLEFILE, as.is = T)
sampleassignments <- strsplit(x = reltable$tumor_wgs_aliquot_id, split = ",")
sampleassignmentsv <- setNames(object = rep(x = reltable$wgs_exclusion_white_gray, lengths(sampleassignments)), nm = unlist(sampleassignments))
#note

resamplingmodel <- read.delim(file = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA_NODRIV.txt", as.is = T)
# resamplingmodel <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Comparison_Permutation-Resampling_simulator_data_noPRAD_CA.txt", as.is = T)
divergent_variants <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary_20210217.txt", as.is = T)


allsums$simnumber_new <- resamplingmodel[match(x = allsums$sampleid, table = resamplingmodel$sampleid), "total_new"]
allsums$ndivergent <- divergent_variants[match(x = allsums$sampleid, table = divergent_variants$sampleid), "nbiallelics"]
allsums$exclusionlist <- sampleassignmentsv[allsums$sampleid]



##### READING AND WRITING OUT THE HITS
# View(allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ]) #ALT
# View(allsums[which(allsums$simratio < 8 & allsums$nparallel > 0 & allsums$snv_slope < 1), ]) #ALT
# View(allsums[which(allsums$nparallel > 0 & allsums$snv_slope < 1 & allsums$is_preferred), ]) #ALT
# View(allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope < 1), ])
# nrow(allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope < 1), ])

write.table(x = allsums, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_20210701.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_20210625.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums[which(allsums$prec > 0 & allsums$rec > 0 & allsums$snv_slope < 1.1), ], file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums[which(allsums$pratio > .01 & allsums$nparallel > 0 & allsums$snv_slope < 1.1), ], file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope <= 1), ],
#             file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210625.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(x = allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope <= 1), ],
            file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210701.txt", quote = F, col.names = T, row.names = F, sep = "\t")

### writing out the specific variants for those samples which passed required filters
parhitsamples <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210701.txt", as.is = T)
# parhitsamples <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210625.txt", as.is = T)
# parhitsamples <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)

# read ACTUAL hits
# parhitsfiles <- paste0("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", parhitsamples$sampleid, "/", parhitsamples$sampleid, "_snv_mnv_infSites_finalhits.txt")
# parhitsfiles <- paste0("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", parhitsamples$sampleid, "/", parhitsamples$sampleid, "_snv_mnv_infSites_finalhits.txt")
parhitsfiles <- paste0("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210702_vafpipeline_out_alphapt1_hetonly/", parhitsamples$sampleid, "/", parhitsamples$sampleid, "_snv_mnv_infSites_finalhits.txt")
parhits <- lapply(X = parhitsfiles, FUN = read.delim, as.is = T)
parhitsdf <- do.call(rbind, parhits)
parhitsdf$sampleid <- rep(parhitsamples$sampleid, sapply(parhits, nrow))


parhitsdf$ref <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$ref))
parhitsdf$alt <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$alt))


# #### 20210625 additional cleaning based on clustering of biallelic parallel hits & flagged mapping/assembly issues 
# #### The germline SV filtering is now included in the early steps and can be omitted here
# dfforclustercheck <- parhitsdf[parhitsdf$sampleid %in% sumtab[sumtab$is_preferred, "samplename"], ]
# tophits <- names(which(table(parhitsdf$sampleid) >= 10))
# # dfforclustercheck <- parhitsdf[parhitsdf$sampleid %in% tophits, c("chr", "start", "ref", "alt", "sampleid")]
# colnames(dfforclustercheck)[2] <- "pos"
# 
# ### run kataegis detection on this stuff
# library(BSgenome.Hsapiens.1000genomes.hs37d5)
# source(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/kataegis/20180309_dpclust3p_kataegis.R")
# source(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/kataegis/20180309_dpclust3p_fastPCF.R")
# 
# ## Rainfall plot
# cumdistv <- setNames(object = c(0, cumsum(as.numeric(seqlengths(BSgenome.Hsapiens.1000genomes.hs37d5)[1:23]))), nm = c(1:22, "X", "Y"))
# dfforclustercheck$cumpos <- dfforclustercheck$pos + cumdistv[dfforclustercheck$chr]
# 
# # do a very lenient kataegis detection, rather overcall than undercall as will be required to overlap with mapping issue region
# identifyKataegis(samplename = "allbiallelics", snvs = dfforclustercheck[!duplicated(dfforclustercheck$cumpos),], outdir = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/", 
#                  maxthresh = 5e3, pstreak = .05, minmutsrange = 4:5,
#                  bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
# katfoci <- read.delim(file = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/allbiallelics_kataegis_annotated.txt", as.is = T)
# 
# 
# dfforclustercheckgr <- GRanges(seqnames = dfforclustercheck$chr, ranges = IRanges(start = dfforclustercheck$pos, width = 1), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
# seqlevelsStyle(dfforclustercheckgr) <- "UCSC"
# 
# 
# # # read mapping/assembly issue annotation and make into GRanges
# # hg19fixtable <- read.delim("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/hg19fix_table.txt.gz", as.is = T)
# # UCSCliftover <- read.delim("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/UCSCliftOverHg38_table.txt.gz", as.is = T)
# # NCBIliftover <- read.delim("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/NCBIReMapHg38_table.txt.gz", as.is = T)
# # hg38ContigDiff <- read.delim("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/hg38ContigDiff.txt.gz", as.is = T)
# # althaplo <- read.delim("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/altHaplotypeAlignment.txt.gz", as.is = T)
# # encodeexclusion <- read.delim("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/ENCFF001TDO.bed.gz", as.is = T, header = F)
# # 
# # hg19fixtablegr <- GRanges(seqnames = hg19fixtable$tName, ranges = IRanges(start = hg19fixtable$tStart, end = hg19fixtable$tEnd))
# # hg38ContigDiffgr <- GRanges(seqnames = hg38ContigDiff$chrom, ranges = IRanges(start = hg38ContigDiff$chromStart, end = hg38ContigDiff$chromEnd))
# # NCBIliftovergr <- GRanges(seqnames = NCBIliftover$tName, ranges = IRanges(start = NCBIliftover$tStart, end = NCBIliftover$tEnd))
# # UCSCliftovergr <- GRanges(seqnames = UCSCliftover$tName, ranges = IRanges(start = UCSCliftover$tStart, end = UCSCliftover$tEnd))
# # althaplogr <- GRanges(seqnames = althaplo$tName, ranges = IRanges(start = althaplo$tStart, end = althaplo$tEnd))
# # encodeexclusiongr <- GRanges(seqnames = encodeexclusion$V1, ranges = IRanges(start = encodeexclusion$V2, end = encodeexclusion$V3))
# # 
# # # annotate df used for kataegis detection (ie. parhits, representative samples only)
# # dfforclustercheck$in_fix <- overlapsAny(query = dfforclustercheckgr, subject = hg19fixtablegr)
# # dfforclustercheck$is_diff <- overlapsAny(query = dfforclustercheckgr, subject = hg38ContigDiffgr)
# # dfforclustercheck$ncbilift_multiple <- countOverlaps(query = dfforclustercheckgr, subject = NCBIliftovergr) > 1
# # dfforclustercheck$ucsclift_multiple <- countOverlaps(query = dfforclustercheckgr, subject = UCSCliftovergr) > 1
# # dfforclustercheck$in_alt <- overlapsAny(query = dfforclustercheckgr, subject = althaplogr)
# # dfforclustercheck$encodeex <- overlapsAny(query = dfforclustercheckgr, subject = encodeexclusiongr)
# # # dfforclustercheck$liftover1kb <- dfforclustercheckgr$hg38clean
# 
# 
# 
# # check common germline CNVs/SVs (notably deletions)
# dbvarcommon1kg <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/dbVar_common_1000g.bed.gz", as.is = T)
# dbvarcommon1kggr <- GRanges(seqnames = dbvarcommon1kg$X.chrom, ranges = IRanges(start = dbvarcommon1kg$chromStart, end = dbvarcommon1kg$chromEnd), type = dbvarcommon1kg$type, freq = as.numeric(gsub(pattern = "ALL_AF=", replacement = "", x = dbvarcommon1kg$frequency)))
# # dbvarcommon1kggr <- dbvarcommon1kggr[!dbvarcommon1kggr$type %in% c("duplication")]
# 
# # dbvarcommonall <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/nstd186.GRCh37.variant_call.tsv.gz", as.is = T, skip = 1)
# # AF <- as.numeric(gsub(pattern = "AF=", replacement = "", x = sapply(X = strsplit(x = dbvarcommonall$allele_frequency, split = ","), FUN = "[", 1)))
# # dbvarcommonallgr <- GRanges(seqnames = dbvarcommonall$chr, ranges = IRanges(start = ifelse(is.na(dbvarcommonall$start), dbvarcommonall$inner_start, dbvarcommonall$start),
# #                                                                             end = ifelse(is.na(dbvarcommonall$stop), dbvarcommonall$inner_stop, dbvarcommonall$stop)), type = dbvarcommonall$variant_call_type, freq = AF, desc = dbvarcommonall$description)
# # dbvarcommonallgr <- dbvarcommonallgr[!dbvarcommonallgr$type %in% c("duplication", "copy number gain")]
# # 
# # dbvardels <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/GRCh37.nr_deletions.common.bed.gz", as.is = T, header = F)
# # dbvarcommonallgr <- GRanges(seqnames = dbvardels$V1, ranges = IRanges(start = dbvardels$V2, end = dbvardels$V3), type = "del")
# # seqlevelsStyle(dbvarcommonallgr) <- "UCSC"
# 
# 
# dfforclustercheck$in_focus <- NA
# dfforclustercheck$in_focus <- katfoci[match(x = paste0(dfforclustercheck$chr, "_", dfforclustercheck$pos), table = paste0(katfoci$chr, "_", katfoci$pos)), "focus"]
# 
# 
# ###NOTE check enrichment = motivation to exclude
# # dfforclustercheck$in_common_sv <- overlapsAny(query = dfforclustercheckgr, subject = dbvarcommon1kggr)
# svoverlaps <- findOverlaps(query = dfforclustercheckgr, subject = dbvarcommon1kggr)
# # svoverlaps <- findOverlaps(query = dfforclustercheckgr, subject = dbvarcommonallgr)
# dfforclustercheck$svtypeoverlap <- "no_common_sv"
# dfforclustercheck[queryHits(svoverlaps), "svtypeoverlap"] <- dbvarcommon1kggr$type[subjectHits(svoverlaps)]
# # dfforclustercheck[queryHits(svoverlaps), "svtypeoverlap"] <- dbvarcommonallgr$type[subjectHits(svoverlaps)]
# table(phasable = dfforclustercheck$is_confirmed, in_sv = dfforclustercheck$svtypeoverlap) # should have 0 in sv as filtered on now
# # table(phasable = dfforclustercheck$is_confirmed, in_sv = !dfforclustercheck$svtypeoverlap %in% "no_common_sv")
# 
# # table(in_focus = !is.na(dfforclustercheck$in_focus), in_sv = dfforclustercheck$in_common_sv)
# # outtab <- table(conf = dfforclustercheck$is_confirmed[dfforclustercheck$is_phaseable], in_sv = dfforclustercheck$in_common_sv[dfforclustercheck$is_phaseable])
# # outtab
# # chisq.test(outtab)
# # table(phasable = dfforclustercheck$is_phaseable, in_sv = dfforclustercheck$in_common_sv)
# # 
# # chisq.test(table(in_focus = !is.na(dfforclustercheck$in_focus), in_sv = dfforclustercheck$in_common_sv))
# # table(dfforclustercheck$is_confirmed[dfforclustercheck$is_phaseable], !is.na(dfforclustercheck$in_focus)[dfforclustercheck$is_phaseable])
# # chisq.test(table(conf = dfforclustercheck$is_confirmed[dfforclustercheck$is_phaseable], clust = !is.na(dfforclustercheck$in_focus)[dfforclustercheck$is_phaseable]))
# # table(conf = dfforclustercheck$is_phaseable, clust = !is.na(dfforclustercheck$in_focus))
# 
# 
# 
# 
# options(scipen = 999)
# ### load all SVs and check
# # allsvfiles <- paste0("/camp/project/proj-emedlab-vanloo/ICGC-structural-variants/", tophits, ".pcawg_consensus_1.6.161116.somatic.sv.vcf.gz")
# # allsvfiles <- allsvfiles[which(file.size(allsvfiles) > 3943)]
# # allsvs <- do.call(rbind, lapply(X = allsvfiles, FUN = function(x) read.delim(x, comment.char = "#", header = F)[,1:2]))
# # allsvs$cumpos <- allsvs$V2 + cumdist[allsvs$V1]
# 
# # dfforclustercheck[duplicated(dfforclustercheck$cumpos) & !is.na(dfforclustercheck$in_focus), ]
# 
# 
# ### sorting stuff
# plotdf <- dfforclustercheck
# plotdf$chr <- factor(x = plotdf$chr, levels = c(1:22, "X"))
# plotdf <- plotdf[order(plotdf$chr, plotdf$pos),]
# 
# plotdf$limd <- log10(do.call(c, by(data = plotdf, INDICES = plotdf$chr, FUN = function(x) c(diff(x$pos), NA))))
# ###
# 
# p1 <- ggplot(data = plotdf, mapping = aes(x = cumpos, y = limd)) +
#   geom_vline(xintercept = cumdistv, color = "grey", size = 0.25) +
#   # geom_point(mapping = aes(colour = svtypeoverlap), show.legend = T, alpha = .75, size = 1)
# geom_point(mapping = aes(colour = is.na(in_focus)), show.legend = F, alpha = .75, size = 1)
# # p1 <- p1 + scale_color_manual(values = c('TRUE' = "grey", 'FALSE' = "red"))#cvect)
# p1 <- p1 + scale_color_brewer(type = "qual", palette = "Accent")
# p1 <- p1 + scale_x_continuous(breaks = c(filter(x = cumdistv, filter = c(0.5, 0.5)))[-24], labels = names(cumdistv)[-24]) + 
#   theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
#   scale_y_continuous(breaks = (0:6), labels = c(1, 10, 100, 1000, 10000, 100000, 1000000), limits = c(0,log10(3e6))) + annotation_logticks(sides = "l")
# p1 <- p1 + labs(y = "log10(inter-mutation distance)", x = "")
# # p1 <- p1 + geom_point(data = allsvs, mapping = aes(x = cumpos, y = 0.1), shape = "I", size = 2, alpha = .01)
# p1
# ggsave(filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Reviewer_biallelic_rainfall_germlineSVs_afterSVfiltering_kafoci.pdf", plot = p1, width = 16, height = 5)
# #   return(NULL)
# # }
# 
# 
# #### cleaning up all the variants in crappy regions (excapt those which are phasing validated)
# katfociranges <- unlist(GRangesList(lapply(X = split(x = katfoci, f = katfoci$focus), FUN = function(x) GRanges(seqnames = x$chr[1], ranges = IRanges(start = min(x$pos), end = max(x$pos)), focus = x$focus[1], nmuts = nrow(x)))))
# parhitsdf$is_clustered <- overlapsAny(query = GRanges(parhitsdf), subject = katfociranges)
# 
# parhitsdf_unfilt <- parhitsdf
# # write.table(x = parhitsdf_unfilt, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210625_wclusteringinfo.txt", quote = F, row.names = F, col.names = T, sep = "\t")
# write.table(x = parhitsdf_unfilt, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210701_wclusteringinfo.txt", quote = F, row.names = F, col.names = T, sep = "\t")
# 
# parhitsdf <- parhitsdf[ifelse(!parhitsdf$is_clustered, T, parhitsdf$is_confirmed), ]
# write.table(x = parhitsdf, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210701_wclusteringinfo_finalset.txt", quote = F, row.names = F, col.names = T, sep = "\t")
# # write.table(x = parhitsdf, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210625_wclusteringinfo_finalset.txt", quote = F, row.names = F, col.names = T, sep = "\t")
# ##### END READING AND WRITING OUT THE HITS
# 
# allsums$nparallel_old <- allsums$nparallel
# allsums$nparallel <- 0
# recountsv <- c(table(parhitsdf$sampleid))
# allsums$nparallel[match(x = names(recountsv), table = allsums$sampleid)] <- recountsv
# 
# # write.table(x = allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope <= 1), ],
# #             file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210625_updatedafterclust.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# write.table(x = allsums[which(allsums$simnumber+allsums$simnumber_divergent >= 0.001 & (allsums$simnumber_new >= 1 | is.na(allsums$simnumber_new)) & allsums$snv_slope <= 1), ],
#             file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_20210701_updatedafterclust.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# 
# 
# ###### end of previous clustering/SV fitlering section - now omitted as included during calling
write.table(x = parhitsdf, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210701_finalset.txt", quote = F, row.names = F, col.names = T, sep = "\t")



allsums <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_20210701.txt", as.is = T)

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(allsums$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(allsums$histology_abbreviation)


# clean set of samples for PLOTTING too
df1 <- allsums[which(allsums$snv_slope <= 1 & allsums$is_preferred), ]
# View(df1[is.na(df1$simnumber_new),])
# df1[is.na(df1$simnumber_new), "simnumber_new"] <- 0
df1[is.na(df1$ndivergent), "ndivergent"] <- 0
df1 <- df1[which(df1$simnumber+df1$simnumber_divergent >= 0.001 & (df1$simnumber_new >= 1 | is.na(df1$simnumber_new))), ]



cv_newsim <- cor.test(df1$simnumber_new, df1$nparallel+df1$ndivergent)
cv_oldsim <- cor.test(df1$simnumber+df1$simnumber_divergent, df1$nparallel+df1$ndivergent)


p1 <- ggplot(data = df1, mapping = aes(x = simnumber+simnumber_divergent, y = nparallel+ndivergent))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_smooth() + 
  geom_point(shape = 21, mapping = aes(fill = histology_abbreviation), color = "black", alpha = .6, size = 5) +
  # geom_point(shape = 21, mapping = aes(fill = exclusionlist), alpha = .6, size = 5) + 
  annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_oldsim$estimate[[1]], digits = 3), collapse = "")) +
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3)) + annotation_logticks() +
  scale_fill_manual(values = cvect, name = "Cancer type") +
  # geom_errorbar(mapping = aes(ymin = lower, ymax = upper)) +
  labs(x = "Average biallelic violations per 1000 simulations (Permutation model)", y = "Number of biallelic variants detected (divergent + parallel)") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4), ylim = c(1, 1e4))
p1
# ggsave(plot = p1, filename = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedPermutation.pdf", width = 15, height = 12)
ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedPermutation.pdf", width = 15, height = 12)




p1 <- ggplot(data = df1, mapping = aes(x = simnumber_new/1000, y = nparallel+ndivergent))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_smooth() + geom_point(shape = 21, mapping = aes(fill = histology_abbreviation), color = "black", alpha = .6, size = 5) + 
  annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_newsim$estimate[[1]], digits = 3), collapse = "")) +
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3)) + annotation_logticks() + scale_fill_manual(values = cvect, name = "Cancer type") + 
  labs(x = "Average biallelic violations per 1000 simulations (Resampling model)", y = "Number of biallelic variants detected (divergent + parallel)") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4), ylim = c(1, 1e4))
p1
ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedResampling_NODRIVERS.pdf", width = 15, height = 12)
write.table(x = df1, file = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulated.txt_NODRIVERS", quote = F, row.names = F, sep = "\t")
# ggsave(plot = p1, filename = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulatedResampling.pdf", width = 15, height = 12)
# write.table(x = df1, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Observed_biallelic_variants_vs_simulated.txt", quote = F, row.names = F, sep = "\t")



### added one to use the estimate based on phasing to the plot
cv_newsim_mod <- cor.test(df1$simnumber_new, ifelse(df1$tot_phaseable>1e4, df1$med, df1$nparallel) + df1$ndivergent)

p1 <- ggplot(data = df1, mapping = aes(x = simnumber_new/1000, y = ifelse(tot_phaseable>1e4, med, nparallel) +ndivergent))
for (idx in 1:8) {
  p1 <- p1 + geom_abline(intercept = c(0,log10(c(0,2,4,8,16,32,64)))[idx], color = paste0("grey", seq(from = 10, to = 90, by = 10))[idx] )
}
p1 <- p1 + geom_smooth() + geom_point(shape = 21, mapping = aes(fill = histology_abbreviation, color = tot_phaseable>1e4, size = tot_hetero), alpha = .6, size = 5) + 
  scale_color_manual(values = c('TRUE' = "red", 'FALSE' = "black"), name = "# phaseable SNVs > 10,000") + scale_size()
  annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_newsim_mod$estimate[[1]], digits = 3), collapse = "")) +
  theme_minimal() + scale_x_log10(labels = 10^(-3:3), breaks = 10^(-3:3)) + scale_y_log10(labels = 10^(0:3), breaks = 10^(0:3)) + annotation_logticks() + scale_fill_manual(values = cvect, name = "Cancer type") + 
  labs(x = "Average biallelic violations per 1000 simulations (Resampling model)", y = "Number of biallelic variants detected/estimated (divergent + parallel)") + coord_fixed(ratio = 1, xlim = c(0.001, 1e4), ylim = c(1, 1e4))
# p1 <- p1 + geom_point(data = df1[df1$tot_phaseable>1e4,], mapping = aes(x = simnumber_new/1000, y = med+ndivergent, fill = histology_abbreviation), alpha = .1)
p1
# ggsave(plot = p1, filename = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/ObservedOrEstimated_biallelic_variants_vs_simulatedResampling.pdf", width = 15, height = 12)
ggsave(plot = p1, filename = "/camp/lab/vanloop/working/demeulj/projects/2016-17_ICGC/infinite_sites/results/figures/ObservedOrEstimated_biallelic_variants_vs_simulatedResampling_NODRIVER.pdf", width = 15, height = 12)




# 
### plot of sim vs obs vs mutation burden
# p1 <- ggplot(data = df1, mapping = aes(x = tot_hetero, y = ifelse(df1$tot_phaseable>1e4, df1$med, df1$nparallel) + df1$ndivergent ))
# p1 <- p1 + geom_point(shape = 21, mapping = aes(fill = histology_abbreviation, color = tot_phaseable>1e4), alpha = .6, size = 5) +
#   scale_color_manual(values = c('TRUE' = "red", 'FALSE' = "black"), name = "# phaseable SNVs > 10,000") +
#   # annotate(geom = "text", x = 0.01, y = 5000, label = paste0("r = ", round(cv_newsim_mod$estimate[[1]], digits = 3), collapse = "")) +
#   theme_minimal() + scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "b") + scale_fill_manual(values = cvect, name = "Cancer type") +
#   labs(x = "Average biallelic violations per 1000 simulations (Resampling model)", y = "Number of biallelic variants detected/estimated (divergent + parallel)")
# # p1 <- p1 + geom_point(data = df1[df1$tot_phaseable>1e4,], mapping = aes(x = simnumber_new/1000, y = med+ndivergent, fill = histology_abbreviation), alpha = .1)
# p1

# 
# cor.test(df1$simnumber_new, ifelse(df1$tot_phaseable>1e4, df1$med, df1$nparallel) + df1$ndivergent)
# cor.test(df1$tot_hetero, ifelse(df1$tot_phaseable>1e4, df1$med, df1$nparallel) + df1$ndivergent)
# cor.test(df1$snvs_total, ifelse(df1$tot_phaseable>1e4, df1$med, df1$nparallel) + df1$ndivergent)
# 
# df1mod <- df1[!is.na(df1$simnumber_new), ]
# df1mod <- df1mod[order(df1mod$tot_testable, decreasing = F), ]
# df1mod$cumsum <- round(cumsum(x = df1mod$tot_hetero) / 1e5, digits = 0)
# 
# v1 <- do.call(rbind, by(data = df1mod, INDICES = df1mod$cumsum, FUN = function(x) c(nvar = x$cumsum[1], sim = sum(x$simnumber_new/1000), obs = sum(ifelse(x$tot_phaseable>1e4, x$med, x$nparallel) + x$ndivergent))))
# 
# plot(v1[,1], log2(v1[,3]/v1[,2]))


# 
# 
# ####### checks comparing 2021 vs 2020 calls
# head(parhitsamples)
# parhits_old <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)
# # parhitsdf_old <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_20210625_wclusteringinfo_finalset.txt", as.is = T)
# inold <- setdiff(x = parhits_old$sampleid, y = parhitsamples$sampleid)
# innew <- setdiff(y = parhits_old$sampleid, x = parhitsamples$sampleid)
# 
# parhits_old[parhits_old$sampleid %in% inold, "nparallel"]
# parhitsamples[parhitsamples$sampleid %in% innew, "nparallel"]
# 
# mergedf <- merge(x = parhits_old, y = parhitsamples, by = "sampleid", all = T, suffixes = c("_old", "_new"))
# 
# p1 <- ggplot(data = mergedf, mapping = aes(x = log10(nparallel_old), y = log10(nparallel_new))) + geom_point() + geom_abline() + annotation_logticks()
# p1
# 
# p1 <- ggplot(data = mergedf, mapping = aes(x = prec_new-prec_old)) + geom_histogram()
# p1
# p1 <- ggplot(data = mergedf, mapping = aes(x = rec_new-rec_old)) + geom_histogram()
# p1
# 
# sum(parhits_old$nparallel)
# sum(parhitsamples$nparallel)
# 
# missingnow <- sumtab[!sumtab$samplename %in% allsums$sampleid, ]
# sort(rowSums(sumtab[!sumtab$samplename %in% allsums$sampleid, c("num_clonal", "num_subclonal")]))
# 
# mergedf[which(mergedf$nparallel_dipl_old - mergedf$nparallel_dipl_new > 100), c("sampleid", "snvs_total", "nparallel_old", "nparallel_new", "prec_old", "prec_new", "rec_old", "rec_new")]

