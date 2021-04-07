#### Figures for InfSites manuscript

# 1b
# load simulations, sum effects of fwd/back, parallel, 3al and plot as points per sample
library(parallel)
library(Biostrings)
library(ggplot2)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

SIMDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown/"
simfiles <- list.files(path = SIMDIR, pattern = "_infsites_backfwd_allelic.txt", full.names = T)
sims <- mclapply(X = simfiles, FUN = read.delim, as.is = T, mc.preschedule = T, mc.cores = 12)

generate_bases_types_trinuc <- function() {
  # generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  types <- c("A", "G", "T", "A", "C", "G")
  types2 <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(trinucleotides, ">", rep(types, rep(16,6))))
  return(list(bases = bases, types = types, types2 = types2, trinucleotides = trinucleotides,
              trinucleotides_mutations = trinucleotides_mutations))
}

generate_isa_violation_df <- function() {
  trinucs <- generate_bases_types_trinuc()
  idxv <- setNames(object = 1:length(trinucs$trinucleotides_mutations), nm = trinucs$trinucleotides_mutations)
  
  mutdf <- expand.grid(mut1 = trinucs$trinucleotides_mutations, mut2 = trinucs$trinucleotides_mutations)
  
  # annotation
  mutdf[, c("reftri1", "alt1")] <- do.call(rbind, strsplit(as.character(mutdf$mut1), split = ">"))
  mutdf[, c("reftri2", "alt2")] <- do.call(rbind, strsplit(as.character(mutdf$mut2), split = ">"))
  mutdf$alttri1 <- DNAStringSet(paste0(substr(mutdf$reftri1, 1,1), mutdf$alt1, substr(mutdf$reftri1, 3,3)))
  mutdf[mutdf$alt1 %in% c("A", "G"), "alttri1"] <- reverseComplement(mutdf[mutdf$alt1 %in% c("A", "G"), "alttri1"])
  mutdf$alttri1 <- as.character(mutdf$alttri1)
  
  mutdf$alttri2 <- DNAStringSet(paste0(substr(mutdf$reftri2, 1,1), mutdf$alt2, substr(mutdf$reftri2, 3,3)))
  mutdf[mutdf$alt2 %in% c("A", "G"), "alttri2"] <- reverseComplement(mutdf[mutdf$alt2 %in% c("A", "G"), "alttri2"])
  mutdf$alttri2 <- as.character(mutdf$alttri2)
  
  mutdf$type <- "empty"
  mutdf$type <- ifelse(mutdf$reftri1 == mutdf$reftri2,
                       ifelse(mutdf$alt1 == mutdf$alt2, "parallel", "third_allele"),
                       ifelse(mutdf$alttri1 == mutdf$reftri2, ifelse(mutdf$alttri2 == mutdf$reftri1, "back", "forward"), mutdf$type))
  mutdf$type <- factor(mutdf$type, levels = c("back", "forward", "parallel", "third_allele", "empty"))
  # mutdf <- mutdf[!(mutdf$type == "third_allele" & as.integer(mutdf$mut1) > as.integer(mutdf$mut2)), ]
  return(mutdf)
}


# read and fix combined histology file
read_histology <- function(histologyfile, melafile = "/srv/shared/vanloo/ICGC_annotations/icgc_melanoma_new_label.txt") {
  histology_all <- read.delim(file = histologyfile, as.is = T)
  melaannots <- read.delim(file = melafile, as.is = T)
  melaannots[melaannots$subtype == "Cutaneous", "subtype"] <- "Cut"
  histology_all$histology_abbreviation <- ifelse(histology_all$histology_abbreviation == "Kidney-RCC",
                                                 ifelse(grepl("papillary", histology_all$histology_tier4), "Kidney-RCC-Pap","Kidney-RCC-Clear"), histology_all$histology_abbreviation)
  histology_all <- merge(x = histology_all, y = melaannots[, c("icgc_aliquot", "subtype")], by.x = "samplename", by.y = "icgc_aliquot", all.x = T)
  histology_all$histology_abbreviation <- ifelse(is.na(histology_all$subtype), histology_all$histology_abbreviation, paste0(histology_all$histology_abbreviation, "-", histology_all$subtype))
  histology_all <- histology_all[, !colnames(histology_all) %in% "subtype"]
  return(histology_all)
}



mutdf <- generate_isa_violation_df()

### has been mod
simssimplif <- as.data.frame(do.call(rbind, lapply(X = sims, FUN = function(x, fact) unlist(by(data = rowSums(x[-9217, c(4,7)]), INDICES = fact, FUN = sum)), fact = mutdf$type)))
simssimplif$sampleid <- gsub(pattern = "_infsites_backfwd_allelic.txt", replacement = "", x = basename(simfiles))
simssimplif <- simssimplif[, grep(pattern = "^empty", x = colnames(simssimplif), invert = T)]

SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/ICGC_annotations/summary_table_combined_annotations_v2.txt"
histology_all <- read_histology(histologyfile = CLEANHISTOLOGYFILE)

simssimplif$mutload <- rowSums(sumtable_whole[match(x = simssimplif$sampleid, table = sumtable_whole$samplename), c("num_clonal", "num_subclonal")])
simssimplif$histology <- sumtable_whole[match(x = simssimplif$sampleid, table = histology_all$samplename), "histology_abbreviation"]

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(simssimplif$histology))), scheme = "tumour.subtype")
names(cvect) <- unique(simssimplif$histology)
cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", '#FF4500', '#FF4500')

simssimpliflog <- simssimplif
simssimpliflog[, -c(25:27)] <- log10(simssimpliflog[, -c(25:27)])
simssimpliflog <- do.call(data.frame,lapply(simssimpliflog, function(x) replace(x, x < 0 ,0)))


p1 <- ggplot(data = simssimpliflog, mapping = aes(x = mutload)) 
# p1 <- p1 + geom_pointrange(mapping = aes(y = parallel.freqal_med, ymin = parallel.freqal_low, ymax = parallel.freqal_hi, color = histology), alpha = .75, show.legend = T)
# p1 <- p1 + geom_pointrange(mapping = aes(y = parallel.freqal_med, ymin = parallel.freqal_low, ymax = parallel.freqal_hi, color = histology), alpha = .75, show.legend = F)
# p1 <- p1 + geom_pointrange(mapping = aes(y = third_allele.freqal_med, ymin = third_allele.freqal_low, ymax = third_allele.freqal_hi, color = histology), alpha = .75, show.legend = F)
# p1 <- p1 + geom_pointrange(mapping = aes(y = forward.freqbf_med, ymin = forward.freqbf_low, ymax = forward.freqbf_hi, color = histology), alpha = .75, show.legend = F)
p1 <- p1 + geom_pointrange(mapping = aes(y = back.freqbf_med, ymin = back.freqbf_low, ymax = back.freqbf_hi, color = histology), alpha = .75, show.legend = F)
p1 <- p1 + scale_x_log10(limits = c(1e4, 3e6)) + annotation_logticks() + scale_color_manual(values = cvect) + scale_y_continuous(breaks = c(0,1,2,3), labels = c(0, 10, 100, 1000), limits = c(0, 3.4))
p1 <- p1 + theme_minimal() + labs(x = "Mutation burden", y = "# back violations")
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190411_figure1b_legend.pdf", plot = p1, width = 6, height = 9)



##### revisit of Fig 1b
library(ggplot2)
library(reshape2)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_all_writefrac/"
i1files <- list.files(path = INDIR, pattern = "_infsites_permut_effgenfrac1.txt", full.names = T, recursive = T)

SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

i1samples <- sub(pattern = "_infsites_permut_effgenfrac1.txt", replacement = "", x = basename(i1files))
i1df <- as.data.frame(do.call(rbind, lapply(X = i1files, FUN = function(x) read.delim(file = x, as.is = T)$counts_tot/1000)))
# i1df$sampleid <- i1samples

reftable <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isavioltypes_reffile.txt", as.is = T)

i1df <- as.data.frame(t(apply(X = i1df, MARGIN = 1, FUN = function(x, f) by(data = x, INDICES = f, FUN = sum), f = reftable$type)))
i1df$sampleid <- i1samples

i1df$totviol <- rowSums(i1df[, 1:4])
i1df$logtotviol <- log10(i1df$totviol)
i1df$backfrac <- i1df$back/i1df$totviol*i1df$logtotviol
i1df$forwarfrac <- i1df$forward/i1df$totviol*i1df$logtotviol
i1df$parfrac <- i1df$parallel/i1df$totviol*i1df$logtotviol
i1df$thirdfrac <- i1df$third_allele/i1df$totviol*i1df$logtotviol

i1df$histology_abbreviation <- sumtable_whole[match(x = i1df$sampleid, table = sumtable_whole$samplename), "histology_abbreviation"]
i1df$mutload <- rowSums(sumtable_whole[match(x = i1df$sampleid, table = sumtable_whole$samplename), c("num_clonal", "num_subclonal")])
i1dfsub <- i1df[i1df$totviol >= 1, ]
i1dfsub$sampleid <- factor(i1dfsub$sampleid, levels = i1dfsub[order(i1dfsub$mutload, decreasing = T), "sampleid"])

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(i1dfsub$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(i1dfsub$histology_abbreviation)

plotdf <- melt(data = i1dfsub[, c("sampleid", "histology_abbreviation", "backfrac", "forwarfrac", "parfrac", "thirdfrac")], id.vars = c("sampleid", "histology_abbreviation"))
plotdf$variable <- factor(plotdf$variable, levels = c("parfrac", "thirdfrac", "forwarfrac", "backfrac"), labels = c("parallel", "third_allele", "forward", "back"))

p1 <- ggplot(data = plotdf, mapping = aes(x = sampleid, y = value, fill = variable)) + geom_col()
p1 <- p1 + geom_segment(mapping = aes(x = as.integer(sampleid)-.5, xend = as.integer(sampleid)+.5, y = -0.02, yend = -0.02, colour = histology_abbreviation), size = 1)
p1 <- p1 + geom_text(mapping = aes(label = ifelse(sampleid %in% c("b07bad52-d44c-4b27-900a-960985bfadec", 
                                                                  "760881cc-c623-11e3-bf01-24c6515278c0",
                                                                  "c9f91ded-3b04-4cd1-8ea6-bbc635a8a4f0",
                                                                  "bd3e88b3-b37c-4641-85fa-d8125ba324ca") & variable == "parallel", "*", ""), y = -0.1))
p1 <- p1 + annotation_logticks(sides = "l") + scale_y_continuous(labels = scales::math_format(10^.x)) + scale_color_manual(values = cvect) + scale_fill_brewer(type = "qual", palette = "Pastel1") + theme_minimal()
p1 <- p1 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_blank()) + labs(y = "# violations expected under uniform",x = "Samples")
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/20190611_figure1b_new.pdf", plot = p1, width = 15, height = 6)


# write.table(x = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/bf_pt_rate_estimates.txt", quote = F, row.names = F, sep = "\t")



# figure 1c
# take 93ff example 
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")

library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(parallel)

genomedf <- as.data.frame(seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
genomedf$offset <- c(0, cumsum(as.numeric(genomedf$seqlengths))[-nrow(genomedf)])

# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "dd67dec6-35dd-4efe-b913-ed4884855365"
sampleid <- "93ff786e-0165-4b02-8d27-806d422e93fc"
sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
sampleid <- "14c5b81d-da49-4db1-9834-77711c2b1d38"

SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

# for (sampleid in sumtable_whole[order(sumtable_whole$num_clonal + sumtable_whole$num_subclonal, decreasing = T), "samplename"]) {

# mclapply(X = sumtable_whole[sumtable_whole$num_clonal + sumtable_whole$num_subclonal > 1e5, "samplename"], FUN = plot_bafpipeline_results, mc.cores = 12, mc.preschedule = T, sumtable_whole = sumtable_whole)  


# plot_bafpipeline_results <- function(sampleid, sumtable_whole) {

rho <- sumtable_whole[sumtable_whole$samplename == sampleid, "purity"]
psit <- sumtable_whole[sumtable_whole$samplename == sampleid, "ploidy"]
psi <- 2*(1-rho)+psit*rho

snvfile <- paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
segmentedvaffile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_baf_logr_summarised_segments.txt"))
vafhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_snv_mnv_infSites_finalhits.txt"))

if (!all(file.exists(snvfile, segmentedvaffile, vafhitsfile))) next

snvs <- readVcf(file = snvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
segmentedvaf <- GRanges(read.delim(file = segmentedvaffile, as.is = T))
snvs <- snvs[!(is.na(info(snvs)$t_alt_count) | is.na(info(snvs)$t_alt_count)) ]

hitsobj <- findOverlaps(query = snvs, subject = segmentedvaf)
any(duplicated(queryHits(hitsobj)))

snvsdf <- data.frame(chr = seqnames(snvs), pos = start(snvs), ref = info(snvs)$t_ref_count, alt = info(snvs)$t_alt_count)[queryHits(hitsobj), ]
# snvsdf[, c("lower", "upper")] <- t(mapply(FUN = HDIofICDF,
#                                           shape1 = snvsdf$alt + 1,
#                                           shape2 = snvsdf$ref + 1,
#                                           MoreArgs = list(ICDFname = qbeta, credMass = .95)))
snvsdf[, c("lower", "upper")] <- t(mcmapply(FUN = HDIofICDF,
                                          shape1 = snvsdf$alt + 1,
                                          shape2 = snvsdf$ref + 1,
                                          MoreArgs = list(ICDFname = qbeta, credMass = .95), mc.preschedule = T, mc.cores = 12))
snvsdf$vaf <- snvsdf$alt/(snvsdf$ref+snvsdf$alt)

snvsdf[, c("mcnlo", "mcn", "mcnhi")] <- snvsdf[, c("lower", "vaf", "upper")]/rho*psi*2^mcols(segmentedvaf)[subjectHits(hitsobj), "mu_logr"]
snvsdf$plotpos <- genomedf[snvsdf$chr, "offset"] + snvsdf$pos


conscn <- read.delim(file = paste0("/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/", sampleid, ".consensus.20170119.somatic.cna.annotated.txt"), as.is = T)
conscn[, c("startoffset", "endoffset")] <- genomedf[conscn$chromosome, "offset"] + conscn[, c("start", "end")]

# read baf pipeline results
vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
# vafhitsdf <- vafhitsdf[vafhitsdf$type == "snv" & vafhitsdf$padj < .05 & vafhitsdf$hetsnpneighbours > 0 & vafhitsdf$nsnphitsoutofrange/vafhitsdf$hetsnpneighbours < .1 & vafhitsdf$nlogrhits/vafhitsdf$nlogrneighbours < .1,  ]

segmentedvaf$mcn <- (segmentedvaf$mu_baf*psi*2^segmentedvaf$mu_logr-(1-rho))/rho
segmentedvaf$startoffset <- genomedf[as.character(seqnames(segmentedvaf)), "offset"] + start(segmentedvaf)
segmentedvaf$endoffset <- genomedf[as.character(seqnames(segmentedvaf)), "offset"] + end(segmentedvaf)
segmentedvafdf <- as.data.frame(segmentedvaf)

snvsdf$outlier <- paste0(snvsdf$chr, "_", snvsdf$pos) %in% paste0(vafhitsdf$chr, "_", vafhitsdf$start)
snvsdf <- snvsdf[order(snvsdf$plotpos, decreasing = F),]



### identify example to show phasing
# phasinghitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20180529_hetSNPs+phasing_out/", sampleid, paste0(sampleid, "_paired_snv-snp_phased_converted.txt"))
# phasinghits <- read.delim(file = phasinghitsfile, as.is = T)

# phasidxs1 <- paste0(phasinghits$chr, "_", phasinghits$pos1)
# phasidxs2 <- paste0(phasinghits$chr, "_", phasinghits$pos2)

# bafidxs <- paste0(snvsdf$chr, "_", snvsdf$pos)

# confirmedhits <- which(phasidxs1 %in% bafidxs | phasidxs2 %in% bafidxs)
snvsdf$is_phased <- paste0(snvsdf$chr, "_", snvsdf$pos) %in% paste0(vafhitsdf$chr, "_", vafhitsdf$start)[vafhitsdf$is_confirmed]
# phasinghits$is_phased <- phasinghits$block %in% phasinghits[confirmedhits, "block"]


chroms <- as.character(1:22)

p2 <- ggplot(data = snvsdf[!snvsdf$outlier & snvsdf$chr %in% chroms & runif(n = nrow(snvsdf)) < 1, ], mapping = aes(x = plotpos, y = mcn))
# p2 <- p2 + geom_point(colour = "grey90", shape = ".", alpha = .75, stroke = 0)
p2 <- p2 + geom_hex(binwidth = c(3.7e6, .05), alpha = .25) + scale_fill_viridis_c(trans = "log10")
# p2 <- p2 + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)
p2 <- p2 + geom_segment(data = conscn[conscn$chromosome %in% chroms,], mapping= aes(x = startoffset, xend = endoffset, y = total_cn, yend = total_cn), colour = "#e34a33", size = 2.5,alpha = 1, position = position_nudge(y = 0.04))
p2 <- p2 + geom_segment(data = conscn[conscn$chromosome %in% chroms,], mapping= aes(x = startoffset, xend = endoffset, y = major_cn, yend = major_cn), colour = "#fdbb84", size = 2.5, alpha = 1, position = position_nudge(y = -0.04))
p2 <- p2 + geom_pointrange(data = snvsdf[snvsdf$outlier & !snvsdf$is_phased & snvsdf$chr %in% chroms, ], mapping = aes(ymin = mcnlo, ymax = mcnhi, y = mcn), colour = "#6baed6", alpha = .9, stroke = 1)
p2 <- p2 + geom_pointrange(data = snvsdf[snvsdf$outlier & snvsdf$is_phased & snvsdf$chr %in% chroms, ], mapping = aes(ymin = mcnlo, ymax = mcnhi, y = mcn), colour = "#08519c", alpha = .9, stroke = 1)
# p2 <- p2 + geom_segment(data = segmentedvafdf, mapping= aes(x = startoffset, xend = endoffset, y = mcn, yend = mcn), colour = "black")
p2 <- p2 + geom_vline(xintercept = genomedf[2:25, "offset"], linetype = "dashed", alpha = .75)
p2 <- p2 + scale_x_continuous(breaks = (genomedf$offset[1:22]+genomedf$offset[2:23])/2, labels = as.character(1:22)) + labs(y = "Copy number", x = "Genomic position")
# p2 <- p2 + scale_x_continuous(breaks = (genomedf$offset[1:5]+genomedf$offset[2:6])/2, labels = as.character(1:5), limits = c(0, 1062541961)) + labs(y = "Copy number", x = "Genomic position")
# p2 <- p2 + scale_fill_distiller(palette= "Spectral", direction=-1) 
p2 <- p2 + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+ylim(c(0,5))#+xlim(c(1e9,2e9)) # xlim(c(1.5e9,2e9))
p2
# & sample(x = rep(c(T,F), times = c(1e5,nrow(snvsdf)-1e5)))
# ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", sampleid, "_bafpipelineplot_hexbin.pdf"), plot = p2, width = 16, height = 5)
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", sampleid, "_bafpipelineplot_hexbin.pdf"), plot = p2, width = 24, height = 5)

# ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/temp/", sampleid, "_bafpipelineplot.png"), plot = p2)
# }




### plot all precision recall numbers
library(ggplot2)
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

# parsum <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", as.is = T)
parsum <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_2021.txt", as.is = T)
sum(parsum$tot_phaseable >= 1e4)
median(parsum[parsum$tot_phaseable >= 1e4, "prec"], na.rm = T)
median(parsum[parsum$tot_phaseable >= 1e4, "rec"], na.rm = T)

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(parsum$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(parsum$histology_abbreviation)

p1 <- ggplot() + geom_point(data = parsum, mapping = aes(x = rec, y = prec, colour = histology_abbreviation, size = log10(nparallel)), position = position_jitter(width = .01, height = .01), alpha = .4)
p1 <- p1 + theme_minimal() + coord_equal(xlim = c(0,1), ylim = c(0,1)) + scale_color_manual(values = cvect) + labs(x = "Recall", y = "Precision")
p1 <- p1 + geom_text(data = data.frame(rec = weighted.mean(x = parsum$rec, w = parsum$nparallel), prec = weighted.mean(x = parsum$prec, w = parsum$nparallel)), mapping = aes(x = rec, y = prec, label = "X"), color = "red")
p1

ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Parallel_Precision-Recall_plot.pdf"), plot = p1, width = 12, height = 9)



##### checking mutation spectra

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")

snvfile <- paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
snvs <- readVcf(file = snvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
snvsdf <- data.frame(chr = seqnames(snvs), pos = start(snvs), ref = as.character(ref(snvs)), alt = as.character(unlist(alt(snvs))))

temp <- get_wg_trinuc_normalisation_factors(bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
plot_mutationspectrum(mutations = snvsdf, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", normalise = T)

vafhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/", sampleid, paste0(sampleid, "_baflogr_gc_results.txt"))
vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
vafhitsdf <- vafhitsdf[vafhitsdf$type == "snv" & vafhitsdf$padj < .05 & vafhitsdf$hetsnpneighbours > 0 & vafhitsdf$nsnphitsoutofrange/vafhitsdf$hetsnpneighbours < .1 & vafhitsdf$nlogrhits/vafhitsdf$nlogrneighbours < .1,  ]

snvsdfsub <- snvsdf[paste0(snvsdf$chr, "_", snvsdf$pos) %in% paste0(vafhitsdf$seqnames, "_", vafhitsdf$start), ]
plot_mutationspectrum(mutations = snvsdfsub, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "ISAviol", normalise = T)

# plot_mutationspectrum(mutations = snvsdf, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "squared", normalise = F)



library(reshape2)

######### final cosine similarity plots & Neff (obs/sim parallel in diploid)
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

parsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)
INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_1plus1_only_writefrac/"
i1files <- list.files(path = INDIR, pattern = "_infsites_totals_effgenfrac0.1.txt", full.names = T, recursive = T)

i1samples <- sub(pattern = "_infsites_totals_effgenfrac0.1.txt", replacement = "", x = basename(i1files))
i1df <- do.call(rbind, lapply(X = i1files, FUN = read.delim, as.is = T))
i1df$sampleid <- rep(i1samples, rep(4, length(i1samples)))
i1df <- reshape(data = i1df, direction = "wide", idvar = "sampleid", timevar = "type")

parsamplesdfm <- merge(x = parsamplesdf, y = i1df, by = "sampleid", all.x = T, all.y = F)

library(ggplot2)

plotdf <- parsamplesdfm[order(parsamplesdfm$med_diploid, decreasing = T), ]
plotdf$sampleid <- factor(x = plotdf$sampleid, levels = plotdf$sampleid)
plotdf <- plotdf[which(plotdf$tot_diploid/plotdf$tot_testable > .10 & plotdf$freq_med.parallel >= 10), ]
plotdf$neff <- plotdf$med_diploid/(plotdf$freq_med.parallel/10)

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(plotdf$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(plotdf$histology_abbreviation)
# cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", '#FF4500', '#FF4500')

exp(weighted.mean(x = log(plotdf$neff), w = ((plotdf$upper_diploid-plotdf$lower_diploid)/(plotdf$freq_med.parallel/10))^-1, na.rm = T))

p1 <- ggplot(data = plotdf, mapping = aes(x = sampleid, y = med_diploid/(freq_med.parallel/10), colour = histology_abbreviation)) + 
  geom_pointrange(mapping = aes(ymax = upper_diploid/(freq_med.parallel/10), ymin = lower_diploid/(freq_med.parallel/10)), show.legend = F) + scale_y_log10() + annotation_logticks(sides = "l")
p1 <- p1 + geom_hline(yintercept = exp(weighted.mean(x = log(plotdf$neff), w = ((plotdf$upper_diploid-plotdf$lower_diploid)/(plotdf$freq_med.parallel/10))^-1, na.rm = T)), alpha =.5, color = "red", linetype = "dashed")
p1 <- p1 + scale_color_manual(values = cvect)
p1 <- p1 + theme_minimal() + labs(y = "# parallel violations estimated / simulated", x = "samples") + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank())
p1

ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/nparobs_over_sim_Neff.pdf"), plot = p1, width = 8, height = 5, useDingbats=FALSE)



####### add in cosine similarities parallel
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

# parsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)
# INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/"
parsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_2021.txt", as.is = T)
INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/"
i1files <- list.files(path = INDIR, pattern = "_parallel_thirdallele_cosinesims_het.txt", full.names = T, recursive = T)

i1samples <- sub(pattern = "_parallel_thirdallele_cosinesims_het.txt", replacement = "", x = basename(i1files))
i1df <- do.call(rbind, lapply(X = i1files, FUN = read.delim, as.is = T))
i1df$sampleid <- rep(i1samples, rep(4, length(i1samples)))
i1df <- reshape(data = i1df[i1df$type == "parallel", ], direction = "wide", idvar = "sampleid", timevar = "comparison")

parsamplesdfm <- merge(x = parsamplesdf, y = i1df, all.x = T, all.y = F)
parsamplesdfm$is_signif <- parsamplesdfm$cosine_hi.all < parsamplesdfm$cosine_low.simulated | 
  parsamplesdfm$cosine_hi.simulated < parsamplesdfm$cosine_low.all

wilcox.test(x = parsamplesdfm[parsamplesdfm$nparallel >= 10, "cosine_med.all"], y = parsamplesdfm[parsamplesdfm$nparallel >= 10, "cosine_med.simulated"], alternative = "less")
quantile(parsamplesdfm[parsamplesdfm$nparallel >= 10, "cosine_med.simulated"], na.rm = T)
quantile(parsamplesdfm[parsamplesdfm$nparallel >= 10, "cosine_med.all"], na.rm = T)
outv <- sapply(X = 1:50, FUN = function(x) wilcox.test(x = parsamplesdfm[parsamplesdfm$nparallel >= x, "cosine_med.all"], y = parsamplesdfm[parsamplesdfm$nparallel >= x, "cosine_med.simulated"], alternative = "less")$p.value)
p1 <- ggplot(data = parsamplesdfm[parsamplesdfm$nparallel >= 1, ]) + geom_density(mapping = aes(x = cosine_med.all), color = "red") + geom_density(mapping = aes(x = cosine_med.simulated), color = "blue")
p1


parsamplesdfm <- parsamplesdfm[which(parsamplesdfm$nparallel >= 10 & parsamplesdfm$is_signif), ]
parsamplesdfm$sampleid <- factor(parsamplesdfm$sampleid, levels = parsamplesdfm$sampleid[order(parsamplesdfm$nparallel, decreasing = T)])

sum(parsamplesdfm$cosine_hi.all < parsamplesdfm$cosine_low.simulated, na.rm = T)
sum(parsamplesdfm$cosine_low.all > parsamplesdfm$cosine_hi.simulated, na.rm = T)

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(parsamplesdfm$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(parsamplesdfm$histology_abbreviation)

p2 <- ggplot(data = parsamplesdfm[parsamplesdfm$is_signif, ], mapping = aes(x = as.integer(sampleid), y = cosine_med.simulated, color = histology_abbreviation)) + 
  geom_pointrange(mapping = aes(ymin = cosine_low.simulated, ymax = cosine_hi.simulated), alpha = .8, show.legend = F)
p2 <- p2 + geom_pointrange(mapping = aes(x = as.integer(sampleid), y = cosine_med.all, ymin = cosine_low.all, ymax = cosine_hi.all), alpha = .2, show.legend = F, shape = 17)
p2 <- p2 + geom_text(mapping = aes(label = ifelse(sampleid == "b07bad52-d44c-4b27-900a-960985bfadec", "*", ""), y = 0.91), color = "black")
p2 <- p2 + scale_color_manual(values = cvect) #+ coord_cartesian(clip = "on", ylim = c(0,1)) #+ scale_alpha_manual(values= c('TRUE' = .8, 'FALSE' = .2))
p2 <- p2 + theme_minimal() + labs(y = "mutation spectrum cosine similarity", x = "samples") + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())
p2

write.table(x = parsamplesdfm, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/parallel_spectra_cosine_sims_new_2021_table.txt", quote = F, row.names = F, sep = "\t")
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/parallel_spectra_cosine_sims_new_2021.pdf"), plot = p2, width = 5, height = 3, useDingbats=FALSE)


####### add in cosine similarities third allele
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

thirdsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary_20210217.txt", as.is = T, sep = "\t")
thirdsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary.txt", as.is = T, sep = "\t")
# parsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", as.is = T)
INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly//"
# INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/"
i1files <- list.files(path = INDIR, pattern = "_parallel_thirdallele_cosinesims_het.txt", full.names = T, recursive = T)

i1samples <- sub(pattern = "_parallel_thirdallele_cosinesims_het.txt", replacement = "", x = basename(i1files))
i1df <- do.call(rbind, lapply(X = i1files, FUN = read.delim, as.is = T))
i1df$sampleid <- rep(i1samples, rep(4, length(i1samples)))
i1df <- reshape(data = i1df[i1df$type == "third_allele", ], direction = "wide", idvar = "sampleid", timevar = "comparison")

# SIMDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_1plus1_only_writefrac/"
SIMDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_het_only_writefrac/"
statsfiles <- file.path(SIMDIR, i1samples, paste0(i1samples, "_generalstats_effgenfrac1.txt"))
# file.exists(statsfiles)
i1df$frac_snvs <- sapply(X = statsfiles, FUN = function(x) if (file.exists(x)) read.delim(file = x, as.is = T, header = F)[4,2] else NA)

thirdsamplesdfm <- merge(x = thirdsamplesdf, y = i1df, all.x = T, all.y = F)
# thirdsamplesdfm <- thirdsamplesdfm[which(thirdsamplesdfm$nbiallelics >= 5 & thirdsamplesdfm$frac_snvs > .01), ] # when using 1plus1 only

wilcox.test(x = thirdsamplesdfm[thirdsamplesdfm$nbiallelics >= 10, "cosine_med.all"], y = thirdsamplesdfm[thirdsamplesdfm$nbiallelics >= 10, "cosine_med.simulated"], alternative = "less")
quantile(thirdsamplesdfm[thirdsamplesdfm$nbiallelics >= 10, "cosine_med.simulated"], na.rm = T)
quantile(thirdsamplesdfm[thirdsamplesdfm$nbiallelics >= 10, "cosine_med.all"], na.rm = T)
p1 <- ggplot(data = thirdsamplesdfm[thirdsamplesdfm$nbiallelics >= 10, ]) + geom_density(mapping = aes(x = cosine_med.all), color = "red") + geom_density(mapping = aes(x = cosine_med.simulated), color = "blue")
p1

thirdsamplesdfm$is_signif <- thirdsamplesdfm$cosine_hi.all < thirdsamplesdfm$cosine_low.simulated | 
  thirdsamplesdfm$cosine_hi.simulated < thirdsamplesdfm$cosine_low.all

thirdsamplesdfm <- thirdsamplesdfm[which(thirdsamplesdfm$nbiallelics >= 10 & thirdsamplesdfm$is_signif), ]
thirdsamplesdfm$sampleid <- factor(thirdsamplesdfm$sampleid, levels = thirdsamplesdfm$sampleid[order(thirdsamplesdfm$nbiallelics, decreasing = T)])

sum(thirdsamplesdfm$cosine_hi.all < thirdsamplesdfm$cosine_low.simulated, na.rm = T)
sum(thirdsamplesdfm$cosine_low.all > thirdsamplesdfm$cosine_hi.simulated, na.rm = T)

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(thirdsamplesdfm$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(thirdsamplesdfm$histology_abbreviation)

p2 <- ggplot(data = thirdsamplesdfm, mapping = aes(x = as.integer(sampleid), y = cosine_med.simulated, color = histology_abbreviation)) +
  geom_pointrange(mapping = aes(ymin = cosine_low.simulated, ymax = cosine_hi.simulated), alpha = .8, show.legend = F)
p2 <- p2 + geom_pointrange(mapping = aes(x = as.integer(sampleid), y = cosine_med.all, ymin = cosine_low.all, ymax = cosine_hi.all), alpha = .2, show.legend = F, shape = 17)
p2 <- p2 + geom_text(mapping = aes(label = ifelse(sampleid == "c9f91ded-3b04-4cd1-8ea6-bbc635a8a4f0", "*", ""), y = 0.35), color = "black")
p2 <- p2 + scale_color_manual(values = cvect) #+ coord_cartesian(clip = "on", ylim = c(0,1)) + scale_alpha_manual(values= c('TRUE' = .8, 'FALSE' = .2))
p2 <- p2 + theme_minimal() + labs(y = "mutation spectrum cosine similarity", x = "samples") + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())
p2

write.table(x = thirdsamplesdfm, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/thrid_allele_spectra_cosine_sims_het_new_2021_table.txt", quote = F, row.names = F, sep = "\t")
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/thrid_allele_spectra_cosine_sims_het_new_2021.pdf"), plot = p2, width = 5, height = 3, useDingbats=FALSE)
# ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/thrid_allele_spectra_cosine_sims_1plus1_subsetCALLS.pdf"), plot = p2, width = 12, height = 5)





##### plot showing all findings of parallel and third allele variant counts + comparison to what is expected
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

thirdsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary_20210217.txt", as.is = T, sep = "\t")
thirdvardf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_20210217.txt", as.is = T)

parsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_2021.txt", as.is = T)
parvardf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_2021.txt", as.is = T)

SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
sumtable_whole$snvburden <- sumtable_whole$num_clonal + sumtable_whole$num_subclonal
sumtable_whole$samplename <- factor(sumtable_whole$samplename, levels = sumtable_whole$samplename[order(sumtable_whole$snvburden, decreasing = T)])

thirdsamplesdf$sampleid <- factor(thirdsamplesdf$sampleid, levels = levels(sumtable_whole$samplename))
parsamplesdf$sampleid <- factor(parsamplesdf$sampleid, levels = levels(sumtable_whole$samplename))

plotdf <- merge(x = parsamplesdf, y = thirdsamplesdf, all.x = T)
plotdf <- plotdf[plotdf$nparallel < plotdf$npar_phased + plotdf$npar_vaf, ]
plotdf$snvburden <- sumtable_whole$snvburden[match(x = plotdf$sampleid, table = sumtable_whole$samplename)]

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(plotdf$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(plotdf$histology_abbreviation)
# cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", '#FF4500', '#FF4500')

plotdf$nbiallelics[which(plotdf$nbiallelics == 0)] <- .7
plotdf$lower[which(plotdf$lower == 0)] <- .7

p1 <- ggplot(data = plotdf, mapping = aes(x = sampleid)) 
p1 <- p1 + geom_col(mapping = aes(y = snvburden, fill = histology_abbreviation), alpha = .5)
p1 <- p1 + geom_linerange(mapping = aes(ymin = lower,ymax = upper), color = "#e41a1c", alpha = .75)
# p1 <- p1 + geom_point(mapping = aes(y = snvburden), color = "grey") 
p1 <- p1 + geom_point(mapping = aes(y = nparallel), color = "#e41a1c", alpha = .5)
p1 <- p1 + geom_point(mapping = aes(y = nbiallelics), color = "#377eb8", alpha = .9)
p1 <- p1 + scale_y_log10() + annotation_logticks(sides = "l") + theme_minimal() + scale_fill_manual(values = cvect)# + coord_cartesian(ylim = c(.7,2.5e6))
p1 <- p1 + theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + labs(x = "Samples", y = "SNV burden")
p1

ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/resultssummary_vaf_phasing_overlap_wide_2021.pdf"), plot = p1, width = 19, height = 3.5, useDingbats=FALSE)


# plot of spheres showing log2-fold over expected == REDONE ######### final cosine similarity plots & Neff (obs/sim parallel in diploid)

library(reshape2)
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

# parsamplesdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)
INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_1plus1_only_writefrac/"

i1df <- do.call(rbind, lapply(X = as.character(plotdf$sampleid), FUN = function(x, indir) {
  ifile <- paste0(indir, x, "/", x, "_infsites_permut_effgenfrac0.1.txt")
  if (file.exists(ifile)) {
    # browser()
    i1 <- read.delim(ifile, as.is = T)
    i1 <- data.frame(t(c(by(data = i1$counts_tot, INDICES = i1$type, FUN = sum))/1000), sampleid = x)
    # i1$sampleid <- x
    } else return(NULL)
  }, indir = INDIR))

thirdvardf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T)

plotdf <- merge(x = parsamplesdf, y = thirdsamplesdf, all.x = T)
plotdf <- plotdf[plotdf$nparallel < plotdf$npar_phased + plotdf$npar_vaf, ]
plotdf$snvburden <- sumtable_whole$snvburden[match(x = plotdf$sampleid, table = sumtable_whole$samplename)]

plotdf2 <- merge(x = plotdf, y = i1df, all.x = T, by = "sampleid")
# plotdf2 <- merge(x = plotdf2, y = i1df, all.x = T, by = "sampleid")
plotdf2$nbiallelics_dipl <- c(by(data = thirdvardf, INDICES = thirdvardf$sampleid, FUN = function(x) sum(x$major_cn == 1 & x$minor_cn == 1, na.rm = T)))[as.character(plotdf2$sampleid)]

plotdf2$parfc <- 10*plotdf2$med_diploid/plotdf2$parallel
plotdf2$thfc <- 10*plotdf2$nbiallelics_dipl/plotdf2$third_allele

# plotdf2[, c("parfc_hi", "parfc_med", "parfc_low")] <- 10*plotdf2$med_diploid/plotdf2[, c("freq_low.x", "freq_med.x", "freq_hi.x")]
# plotdf2[, c("thfc_hi", "thfc_med", "thfc_low")] <- 10*plotdf2$nbiallelics_dipl/plotdf2[, c("freq_low.y", "freq_med.y", "freq_hi.y")]

p2 <- ggplot(data = plotdf2, mapping = aes(x = sampleid)) + geom_point(mapping = aes(y = .1, size = log2(thfc), alpha = log10(nbiallelics_dipl)), color = "#377eb8")
p2 <- p2 + geom_point(mapping = aes(y = -.1, size = log2(parfc), alpha = log10(npar_phased)), color = "#e41a1c")
p2 <-  p2 + scale_size_continuous(limits = c(1,6)) + coord_cartesian(ylim = c(-.15, .15))
p2 <- p2 + theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank())
p2
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/resultssummary_vaf_phasing_overlap_obsoverexp.pdf"), plot = p2, width = 18, height = 3.6, useDingbats=FALSE)


# get the numbers
i1df <- do.call(rbind, lapply(X = union(thirdsamplesdf$sampleid, parsamplesdf$sampleid), FUN = function(x, indir) {
  ifile <- paste0(indir, x, "/", x, "_infsites_permut_effgenfrac0.1.txt")
  if (file.exists(ifile)) {
    # browser()
    i1 <- read.delim(ifile, as.is = T)
    i1 <- data.frame(t(c(by(data = i1$counts_tot, INDICES = i1$type, FUN = sum))/1000), sampleid = x)
    # i1$sampleid <- x
  } else return(NULL)
}, indir = INDIR))

i1df$nbial_dipl <- c(by(data = thirdvardf, INDICES = thirdvardf$sampleid, FUN = function(x) sum(x$major_cn == 1 & x$minor_cn == 1, na.rm = T)))[as.character(i1df$sampleid)]
i1df$npar_dipl <- parsamplesdf[match(x = i1df$sampleid, table = parsamplesdf$sampleid), "med_diploid"]

i1df$parfc <- 10*i1df$npar_dipl/i1df$parallel
i1df$thfc <- 10*i1df$nbial_dipl/i1df$third_allele

## 
median(i1df$thfc[which(i1df$third_allele > 1)], na.rm = T)
median(i1df$parfc[which(i1df$parallel > 1)], na.rm = T)




# p2 <- ggplot(data = plotdf2, mapping = aes(x = sampleid)) +   geom_point(mapping = aes(y = .1, size = log2(thfc_med), alpha = log10(nbiallelics_dipl)), color = "#377eb8")
# p2 <- p2 + geom_point(mapping = aes(y = -.1, size = log2(parfc_med), alpha = log10(npar_phased)), color = "#e41a1c")
# p2 <-  p2 + scale_size_continuous(limits = c(1,6)) + coord_cartesian(ylim = c(-.15, .15))
# p2 <- p2 + theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank())
# p2
# ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/resultssummary_vaf_phasing_overlap_obsoverexp.pdf"), plot = p1, width = 12.5, height = 2.5, useDingbats=FALSE)
# 
# # p2 <- p2 + geom_point(mapping = aes(y = -.1, size = log2(parfc_hi)), shape = 1, alpha = .5)
# # p2 <- p2 +geom_point(mapping = aes(y = -.1, size = log2(parfc_low)), shape = 1, alpha = .5)
# 
# p2 <- ggplot(data = plotdf2, mapping = aes(x = sampleid)) + geom_point(mapping = aes(y = thfc_med, alpha = log10(nbiallelics_dipl)), color = "blue") +
#   geom_point(mapping = aes(y = parfc_med, alpha = log10(npar_phased)), color = "red") + scale_y_log10() + annotation_logticks(sides = "l")
# p2
# 
# plotdf2$nbiallelics_dipl[which(plotdf2$nbiallelics_dipl == 0)] <- .1
# 
# p2 <- ggplot(data = plotdf2, mapping = aes(x = sampleid)) + geom_tile(mapping = aes(y = .5, width = 1, height = 1, fill = log2(parfc_med), alpha = log10(npar_phased))) +
#   geom_tile(mapping = aes(y = -.5, width = 1, height = 1, fill = log2(thfc_med), alpha = log10(nbiallelics_dipl)))
# p2 <-  p2 + theme_minimal() + scale_fill_continuous(na.value = 'white') + theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank())
# p2

#### plots showing all parallel/third allele variants in hypermut colorectals
library(BSgenome.Hsapiens.1000genomes.hs37d5)
sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
sampleid <- "14c5b81d-da49-4db1-9834-77711c2b1d38"



get_mutspectrum <- function(mutations, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5) {
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  if (!"trinuc" %in% colnames(mutations))
    mutations$trinuc <- as.character(get_trinuc_context(mutations = mutations, bsgenome = bsgenome))
  
  # reverse complement and data augmentation
  revcomp <- data.frame(alt = as.character(complement(DNAStringSet(mutations$alt))),
                        trinuc = as.character(reverseComplement(DNAStringSet(mutations$trinuc))),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                            paste0(mutations$trinuc, ">", mutations$alt),
                            paste0(revcomp$trinuc, ">", revcomp$alt)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations"]])
  
  # compute frequencies and normalised probabilities
  muttype_freq <- c(table(mut_full))
  return(muttype_freq)
}



get_mutspectrum_biallelic <- function(mutations, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5) {
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  if (!"trinuc" %in% colnames(mutations))
    mutations$trinuc <- as.character(get_trinuc_context(mutations = mutations, bsgenome = bsgenome))
  
  altbases <- split(cbind(mutations$alt1, mutations$alt2), f = 1:nrow(mutations))
  alt <- ifelse(mutations$ref == "C", ifelse(sapply(altbases, setequal, y = c("A","G")), "A+G", 
                                             ifelse(sapply(altbases, setequal, y = c("T","G")), "G+T", "T+A")),
                ifelse(mutations$ref == "G", ifelse(sapply(altbases, setequal, y = c("T","C")), "A+G", 
                                                    ifelse(sapply(altbases, setequal, y = c("A","C")), "G+T", "T+A")),
                       ifelse(mutations$ref == "T", ifelse(sapply(altbases, setequal, y = c("A","C")), "A+C", 
                                                           ifelse(sapply(altbases, setequal, y =c("C","G")), "C+G", "G+A")),
                              ifelse(sapply(altbases, setequal, y = c("T","G")), "A+C", 
                                     ifelse(sapply(altbases, setequal, y =c("G","C")), "C+G", "G+A")))
                ))
  
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"),
                            paste0(mutations$trinuc, ">", alt), 
                            paste0(as.character(reverseComplement(DNAStringSet(mutations$trinuc))), ">", alt)),
                     levels = paste0(base_type_trinuc_info[["trinucleotides_mutations"]], "+", rep(c("G", "T", "A", "C", "G", "A"), rep(16,6))))
  
  # compute frequencies and normalised probabilities
  muttype_freq <- c(table(mut_full))
  return(muttype_freq)
}



# get the trinucleotide context of a set of point mutations
get_trinuc_context <- function(mutations, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, size = 1) {
  sequences <- getSeq(x = bsgenome, names = mutations$chr, start = mutations$start - size, end = mutations$end + size, strand = "+")
  return(sequences)
}


allobshits <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)

allobshits <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt", as.is = T)
sampleobshits <- allobshits[allobshits$sampleid == sampleid, ]

parmutstab <- as.data.frame(get_mutspectrum(mutations = sampleobshits, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5))
colnames(parmutstab) <- "cnt"
parmutstab$mut1 <- factor(rownames(parmutstab), levels = rownames(parmutstab))
parmutstab$type <- paste0(substr(parmutstab$mut1,2,2), substr(parmutstab$mut1,4,5))
parmutstab$trinuc <- substr(parmutstab$mut1,1,3)

bases <- generate_bases_types_trinuc()

# parallel
p1 <- ggplot(data = parmutstab, mapping = aes(x = mut1)) + 
  geom_col(mapping = aes(y = cnt, fill = type), show.legend = F, alpha = .6)
p1 <- p1 + geom_segment(data = data.frame(start = seq(0.65, 96, 16), end = seq(16.35, 97, 16), type = bases$types2),
                        mapping = aes(x = start, xend = end, y = 1.25*max(parmutstab$cnt), yend = 1.25*max(parmutstab$cnt), color = type), size = 3, show.legend = F, alpha = .6)
p1 <- p1 + geom_text(data = data.frame(type = bases$types2, pos = seq(8, 96, 16)),
                     mapping = aes(x = pos, y = 1.3*max(parmutstab$cnt), color = type, label = type),
                     show.legend = F, family = "mono")
p1 <- p1 + theme_bw() + scale_x_discrete(labels = parmutstab$trinuc) + scale_fill_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
  scale_color_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
  coord_cartesian(ylim = c(0, 1.15*max(parmutstab$cnt)), clip="off") + ylab(label = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, family = "mono"), panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(2,.5,.5,.5), "lines"))
p1




p1 <- ggplot(data = parmutstab, mapping = aes(x = fulltype)) + 
  geom_col(mapping = aes(y = cnt, fill = typesplit), position = "stack", show.legend = F, alpha = .6)
p1
p1 <- p1 + geom_segment(data = data.frame(start = seq(0.65, 96, 16), end = seq(16.35, 97, 16), type = paste0(bases$types2, "+", c("G", "T", "A", "C", "G", "A"))),
                        mapping = aes(x = start, xend = end, y = 1.25*max(parmutstab$cnt), yend = 1.25*max(parmutstab$cnt), color = type), size = 3, show.legend = F, alpha = .6)
p1 <- p1 + geom_text(data = data.frame(type = paste0(bases$types2, "+", c("G", "T", "A", "C", "G", "A")), pos = seq(8, 96, 16)),
                     mapping = aes(x = pos, y = 1.3*max(parmutstab$cnt), color = type, label = type),
                     show.legend = F, family = "mono")
p1 <- p1 + theme_bw() + scale_x_discrete(labels = plotdf2$trinuc) + scale_fill_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
  scale_color_manual(values = c("C>A+G" = "#15bcee","C>G+T" = "#000000","C>T+A" = "#e32926","T>A+C" = "#999999","T>C+G" = "#a1ce63","T>G+A" = "#ebc6c4")) +
  coord_cartesian(ylim = c(0, 1.15*max(plotdf2$pobs_hi)), clip="off") + ylab(label = "Probability") +
  theme(axis.text.x = element_text(angle = 90, family = "mono"), panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(2,.5,.5,.5), "lines"))
# p1

outfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_thirdallele_snv_spectra_", simsubset))
ggsave(filename = paste0(outfile, ".pdf"), plot = p1, width = 16, height = 4, useDingbats=FALSE)

