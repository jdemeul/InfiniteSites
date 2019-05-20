### linking observed spectra to predicted spectra (1+1 regions) both parallel & third alleles
library(ggplot2)
library(VariantAnnotation)
library(BSgenome.Hsapiens.1000genomes.hs37d5)


sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"


##### checking mtuation spectra

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")

snvfile <- paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
snvs <- readVcf(file = snvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
snvsdf <- data.frame(chr = seqnames(snvs), pos = start(snvs), ref = as.character(ref(snvs)), alt = as.character(unlist(alt(snvs))))

temp <- get_wg_trinuc_normalisation_factors(bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
plot_mutationspectrum(mutations = snvsdf, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "all2", normalise = T)

vafhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_snv_mnv_infSites_finalhits.txt"))
vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
# vafhitsdf <- vafhitsdf[vafhitsdf$type == "snv" & vafhitsdf$padj < .05 & vafhitsdf$hetsnpneighbours > 0 & vafhitsdf$nsnphitsoutofrange/vafhitsdf$hetsnpneighbours < .1 & vafhitsdf$nlogrhits/vafhitsdf$nlogrneighbours < .1,  ]

snvsdfsub <- vafhitsdf[, c("chr", "start", "ref", "alt")]
colnames(snvsdfsub)[2] <- "pos"
plot_mutationspectrum(mutations = snvsdfsub, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "ISAviol2", normalise = T)

# plot_mutationspectrum(mutations = snvsdf, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "squared", normalise = F)

