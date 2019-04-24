#### reboot on 21/08/2018
#### rewrite the SNP-SNV

### Run heterozygous SNP collection and phasing pipelines

library(readr)
library(GenomicRanges)
library(parallel)
library(ggplot2)
library(MASS)
library(VGAM)
library(ROCR)
library(boot)
library(VariantAnnotation)
# library(rslurm)


# genome
library(BSgenome.Hsapiens.1000genomes.hs37d5)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/getCleanHetSNPs.R")
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/code/20190415_GCbyAF_functions_modforInfSites.R")

RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
REFALLELESDIR <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_loci/"
LOGRDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_logR/"
RHOPSI <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/consensus.20170217.purity.ploidy.txt.gz"
CNDIR <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"


####

# rhopsi
rhopsi <- read.delim(file = RHOPSI, as.is = T)

NCORES <- 8

TEMPDIR <- "/home/jdemeul/temp/"
ALLELECOUNTSDIR <- "/srv/shared/vanloo/ICGC_copynumber/battenberg_raw_files/allelecounts/"
SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/"

BAFDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_rerun_005_input_to_finalConsCopynum_BAFphased/"

# CNDIR <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"
BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out/"
# BAMDIR <- "/srv/shared/vanloo/ICGC/"
PHASINGDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20181021_hetSNPs_all+phasing_out/"



### load one time only
releasetable <- read_pcawg_release_table(release_table_file = RELEASETABLEFILE)
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
# load reference alleles
refalleles_gr <- load_1000G_reference_alleles(refallelesdir = REFALLELESDIR)

# immune loci
immune_loci <- GRanges(seqnames = c(14, 7, 7, 14, 2, 22, 6), ranges = IRanges(start = c(22090057, 141998851, 38279625, 106032614, 89156674, 22380474, 28477797),
                                                                              end = c(23021075, 142510972, 38407656, 107288051, 90274235, 23265085, 33448354)), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))


####################
# start with sample ID (i.e. tumour wgs aliquot ID)
# sampleid <- releasetable[1, "tumor_wgs_aliquot_id"]
# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
sampleid <- "1ac15380-04a2-42dd-8ade-28556a570e80"
# sampleid <- "04aa6b77-8074-480c-872e-a1a47afa5314"
####################


run_baflogr_pipeline <- function(sampleid) {
  
  print(paste0("Loading data for sample ", sampleid))
  
  options(bitmapType = "cairo")
  
  sampledir <- file.path(BASEOUT, sampleid)
  #######
  samplerhopsi <- rhopsi[rhopsi$samplename == sampleid, c("purity", "ploidy"), drop = T]
  inferred_sex <- sumtable_whole[sumtable_whole$samplename == sampleid, "inferred_sex", drop = T]
  
  allelecounts_tum_gr <- get_tum_allecounts_gr(allelecountsdir = ALLELECOUNTSDIR,
                                               sampleid = sampleid,
                                               bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5,
                                               reference_alleles = refalleles_gr,
                                               tempdir = TEMPDIR)

  if (is.null(allelecounts_tum_gr) | !dir.exists(sampledir) )
    return(NULL)

  segments_gr <- get_consensus_cn(sampleid = sampleid,
                                  cndir = CNDIR,
                                  bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
  
  if (inferred_sex == "male")
    segments_gr <- trim_XY_to_PAR(segments_gr, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
  
  phasedbaf_gr <- get_phased_BAF(bafdir = BAFDIR,
                                 sampleid = sampleid,
                                 bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5,
                                 allelecounts = allelecounts_tum_gr)
  
  rm(allelecounts_tum_gr)
  
  gccorrlogr_gr <- get_GCcorr_logr(logrdir = LOGRDIR,
                                   sampleid = sampleid,
                                   # segments_gr = segments_gr, 
                                   bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
  
  snvs_vcf <- get_snv_mnvs(sampleid = sampleid, 
                           snvdir = SNVMNVINDELDIR)
  

  ###### check for load
  if (any(is.null(segments_gr), is.null(phasedbaf_gr), is.null(gccorrlogr_gr), is.null(snvs_vcf)))
    return(NULL)
  ######
  
  segments_gr <- summarise_baflogr_segments(sampleid = sampleid,
                                            sampledir = sampledir, 
                                            phasedbaf = phasedbaf_gr, 
                                            logr = gccorrlogr_gr, 
                                            segments_gr = segments_gr, 
                                            rhopsi = samplerhopsi)
  
  
  # debug(test_clean_sites)
  print(paste0("Testing loci for sample ", sampleid))
  hitvariants <- test_clean_sites(sampledir = sampledir,
                                  sampleid = sampleid,
                                  segments_gr = segments_gr,
                                  phasedbaf_gr = phasedbaf_gr,
                                  logr = gccorrlogr_gr,
                                  snvs_vcf = snvs_vcf,
                                  bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5,
                                  ncores = NCORES, 
                                  presched = T,
                                  subsample_optim = T,
                                  pseudocountrange = c(50, 1000),
                                  recompute_pseudocoverage = F,
                                  immune_loci = immune_loci)
  
  
  print(paste0("Annotating with phasing and deriving violations for sample ", sampleid))
  # debug(get_metrics)
  # debug(get_performance_metrics)
  # debug(call_parallel_violations)
  finhits <- call_parallel_violations(sampleid = sampleid,
                                      sampledir = sampledir,
                                      phasingdir = PHASINGDIR,
                                      nboot = 100,
                                      alpha = .01)
  
  
  # debug(annotate_gc_hits)
  # print(paste0("Annotating loci for sample ", sampleid))
  # hitvariants_annotated <- annotate_gc_hits(hitvariants = hitvariants$hitvariants,
  #                                           pseudocount_calibr = hitvariants$pseudocount_calibr,
  #                                           sampledir = sampledir, 
  #                                           sampleid = sampleid,
  #                                           baf = phasedbaf_gr,
  #                                           logr = gccorrlogr_gr,
  #                                           snvs = snvs_vcf,
  #                                           seg_gr = segments_gr,
  #                                           conversionlength = 1e3,
  #                                           testwindow = 2e4,
  #                                           fdr = .05,
  #                                           plotting = T)
  
  return(NULL)
}




#### rslurm submission command
rslurmdf <- releasetable[, "tumor_wgs_aliquot_id", drop = F]
colnames(rslurmdf) <- "sampleid"

# debug(run_baflogr_pipeline)
# debug(test_clean_sites)
# debug(annotate_gc_hits)
# run_baflogr_pipeline(sampleid = rslurmdf[4,])
run_baflogr_pipeline(sampleid = sampleid)
# baflogrjob <- slurm_apply(f = run_baflogr_pipeline, params = rslurmdf[,, drop = F], jobname = "baflogr_run4", nodes = 463, cpus_per_node = 6, add_objects = ls(),
#                           pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(), submit = T)

# baflogrjob <- slurm_apply(f = run_baflogr_pipeline, params = rslurmdf[,, drop = F], jobname = "baflogr_run4", nodes = 100, cpus_per_node = 6, add_objects = ls(),
#                           pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(), submit = T)
