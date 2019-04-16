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
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(trinucleotides, ">", rep(types, rep(16,6))))
  return(list(bases = bases, types = types, trinucleotides = trinucleotides,
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

simssimplif <- as.data.frame(do.call(rbind, mclapply(X = sims, FUN = function(x, fact) unlist(by(data = x[-9217, 3:8], INDICES = fact, FUN = colSums)), fact = mutdf$type)))
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



# figure 1c
# take deb9f example
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")

library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)

genomedf <- as.data.frame(seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
genomedf$offset <- c(0, cumsum(as.numeric(genomedf$seqlengths))[-nrow(genomedf)])

sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
sampleid <- "dd67dec6-35dd-4efe-b913-ed4884855365"
# sampleid <- "2790b964-63e3-49aa-bf8c-9a00d3448c25"

SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

# for (sampleid in sumtable_whole[order(sumtable_whole$num_clonal + sumtable_whole$num_subclonal, decreasing = T), "samplename"]) {

# mclapply(X = sumtable_whole[sumtable_whole$num_clonal + sumtable_whole$num_subclonal > 1e5, "samplename"], FUN = plot_bafpipeline_results, mc.cores = 12, mc.preschedule = T, sumtable_whole = sumtable_whole)  


# plot_bafpipeline_results <- function(sampleid, sumtable_whole) {

rho <- sumtable_whole[sumtable_whole$samplename == sampleid, "purity"]
psit <- sumtable_whole[sumtable_whole$samplename == sampleid, "ploidy"]
psi <- 2*(1-rho)+psit*rho

snvfile <- paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
segmentedvaffile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/", sampleid, paste0(sampleid, "_baf_logr_summarised_segments.txt"))
vafhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/", sampleid, paste0(sampleid, "_baflogr_gc_results.txt"))

if (!all(file.exists(snvfile, segmentedvaffile, vafhitsfile))) next

snvs <- readVcf(file = snvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
segmentedvaf <- GRanges(read.delim(file = segmentedvaffile, as.is = T))
snvs <- snvs[!(is.na(info(snvs)$t_alt_count) | is.na(info(snvs)$t_alt_count)) ]

hitsobj <- findOverlaps(query = snvs, subject = segmentedvaf)
any(duplicated(queryHits(hitsobj)))

snvsdf <- data.frame(chr = seqnames(snvs), pos = start(snvs), ref = info(snvs)$t_ref_count, alt = info(snvs)$t_alt_count)[queryHits(hitsobj), ]
snvsdf[, c("lower", "upper")] <- t(mapply(FUN = HDIofICDF,
                                          shape1 = snvsdf$alt + 1,
                                          shape2 = snvsdf$ref + 1,
                                          MoreArgs = list(ICDFname = qbeta, credMass = .95)))
snvsdf$vaf <- snvsdf$alt/(snvsdf$ref+snvsdf$alt)

snvsdf[, c("mcnlo", "mcn", "mcnhi")] <- snvsdf[, c("lower", "vaf", "upper")]/rho*psi*2^mcols(segmentedvaf)[subjectHits(hitsobj), "mu_logr"]
snvsdf$plotpos <- genomedf[snvsdf$chr, "offset"] + snvsdf$pos


conscn <- read.delim(file = paste0("/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/", sampleid, ".consensus.20170119.somatic.cna.annotated.txt"), as.is = T)
conscn[, c("startoffset", "endoffset")] <- genomedf[conscn$chromosome, "offset"] + conscn[, c("start", "end")]

# read baf pipeline results
vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
vafhitsdf <- vafhitsdf[vafhitsdf$type == "snv" & vafhitsdf$padj < .05 & vafhitsdf$hetsnpneighbours > 0 & vafhitsdf$nsnphitsoutofrange/vafhitsdf$hetsnpneighbours < .1 & vafhitsdf$nlogrhits/vafhitsdf$nlogrneighbours < .1,  ]

segmentedvaf$mcn <- (segmentedvaf$mu_baf*psi*2^segmentedvaf$mu_logr-(1-rho))/rho
segmentedvaf$startoffset <- genomedf[as.character(seqnames(segmentedvaf)), "offset"] + start(segmentedvaf)
segmentedvaf$endoffset <- genomedf[as.character(seqnames(segmentedvaf)), "offset"] + end(segmentedvaf)
segmentedvafdf <- as.data.frame(segmentedvaf)

snvsdf$outlier <- paste0(snvsdf$chr, "_", snvsdf$pos) %in% paste0(vafhitsdf$seqnames, "_", vafhitsdf$start)
snvsdf <- snvsdf[order(snvsdf$plotpos, decreasing = F),]



### identify example to show phasing
phasinghitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/", sampleid, paste0(sampleid, "_paired_snv-snp_phased_converted.txt"))
phasinghits <- read.delim(file = phasinghitsfile, as.is = T)

phasidxs1 <- paste0(phasinghits$chr, "_", phasinghits$pos1)
phasidxs2 <- paste0(phasinghits$chr, "_", phasinghits$pos2)

bafidxs <- paste0(snvsdf$chr, "_", snvsdf$pos)

# confirmedhits <- which(phasidxs1 %in% bafidxs | phasidxs2 %in% bafidxs)
snvsdf$is_phased <- bafidxs %in% c(phasidxs1, phasidxs2)
# phasinghits$is_phased <- phasinghits$block %in% phasinghits[confirmedhits, "block"]





chroms <- as.character(1:22, "X", "Y")

p2 <- ggplot(data = snvsdf[!snvsdf$outlier & snvsdf$chr %in% chroms & runif(n = nrow(snvsdf)) < 1, ], mapping = aes(x = plotpos, y = mcn))
p2 <- p2 + geom_point(colour = "grey90", shape = 20, alpha = .75, stroke = 0)
# p2 <- p2 + geom_hex(binwidth = c(5e6, .05), alpha = .5)
# p2 <- p2 + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)
p2 <- p2 + geom_segment(data = conscn[conscn$chromosome %in% chroms,], mapping= aes(x = startoffset, xend = endoffset, y = total_cn, yend = total_cn), colour = "#33a02c", size = 2.5,alpha = 1, position = position_nudge(y = 0.04))
p2 <- p2 + geom_segment(data = conscn[conscn$chromosome %in% chroms,], mapping= aes(x = startoffset, xend = endoffset, y = major_cn, yend = major_cn), colour = "#b2df8a", size = 2.5, alpha = 1, position = position_nudge(y = -0.04))
p2 <- p2 + geom_pointrange(data = snvsdf[snvsdf$outlier & !snvsdf$is_phased & snvsdf$chr %in% chroms, ], mapping = aes(ymin = mcnlo, ymax = mcnhi, y = mcn), colour = "#a6cee3", alpha = .9, stroke = 1)
p2 <- p2 + geom_pointrange(data = snvsdf[snvsdf$outlier & snvsdf$is_phased & snvsdf$chr %in% chroms, ], mapping = aes(ymin = mcnlo, ymax = mcnhi, y = mcn), colour = "#1f78b4", alpha = .9, stroke = 1)
# p2 <- p2 + geom_segment(data = segmentedvafdf, mapping= aes(x = startoffset, xend = endoffset, y = mcn, yend = mcn), colour = "black")
p2 <- p2 + geom_vline(xintercept = genomedf[2:25, "offset"], linetype = "dashed", alpha = .75)
p2 <- p2 + scale_x_continuous(breaks = (genomedf$offset[1:5]+genomedf$offset[2:6])/2, labels = as.character(1:5))
# p2 <- p2 + scale_fill_distiller(palette= "Spectral", direction=-1) 
p2 <- p2 + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+ylim(c(0,5))#+xlim(c(1e9,2e9)) # xlim(c(1.5e9,2e9))
p2
# & sample(x = rep(c(T,F), times = c(1e5,nrow(snvsdf)-1e5)))
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", sampleid, "_bafpipelineplot.pdf"), plot = p2, width = 16, height = 5)

ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/temp/", sampleid, "_bafpipelineplot.png"), plot = p2)
# }



##### checking mtuation spectra

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



