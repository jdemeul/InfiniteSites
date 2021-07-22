###!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### precancer ASCAT preprocessing
TUMOURNAME <- as.character(args[1])
# TUMOURNAME <- "S094"
print(TUMOURNAME)


### Run heterozygous SNP collection and phasing pipelines

library(readr)
library(GenomicRanges)
library(rtracklayer)
# library(parallel)
library(ggplot2)
library(MASS)
library(VGAM)
library(ROCR)
library(boot)
library(VariantAnnotation)
# library(rslurm)
library(metap)
# library(rtracklayer)


# genome
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.NCBI.GRCh38)

source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")
source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/getCleanHetSNPs.R")
source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/code/code/20190415_GCbyAF_functions_modforInfSites_rerun20210629.R")
source("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/code/code/check_germline_artefact.R")

RELEASETABLEFILE <- "/camp/project/proj-vanloo-secure/PCAWG/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
# REFALLELESDIR <- "/camp/project/proj-emedlab-vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_loci/"
# LOGRDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_logR/"
REFALLELESDIR <- "/camp/project/proj-emedlab-vanloo/mtarabichi/G1000/GRCh37/Simplified/"
LOGRDIR <- "/camp/project/proj-emedlab-vanloo/mtarabichi/Battenberg/BB_pcawg_fullrerun/"
# LOGRDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/partial_BB_rerun/"
RHOPSI <- "/camp/project/proj-vanloo-secure/PCAWG/ICGC_consensus_copynumber/20170119_release/consensus.20170217.purity.ploidy.txt.gz"
CNDIR <- "/camp/project/proj-vanloo-secure/PCAWG/ICGC_consensus_copynumber/20170119_release/"


hg19tohg38chainfile <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/hg19ToHg38.over.chain"
hg19tohg38chain <- import.chain(con = hg19tohg38chainfile)

# Annotate/remove common germline CNVs/SVs (notably deletions)
dbvarcommon1kg <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/hg19-38fixes/dbVar_common_1000g.bed.gz", as.is = T)
dbvarcommon1kg <- GRanges(seqnames = dbvarcommon1kg$X.chrom, ranges = IRanges(start = dbvarcommon1kg$chromStart, end = dbvarcommon1kg$chromEnd), type = dbvarcommon1kg$type, freq = as.numeric(gsub(pattern = "ALL_AF=", replacement = "", x = dbvarcommon1kg$frequency)))
seqlevelsStyle(dbvarcommon1kg) <- "Ensembl"

####

# rhopsi
rhopsi <- read.delim(file = RHOPSI, as.is = T)

NCORES <- 5

# TEMPDIR <- "/home/jdemeul/temp/"
# ALLELECOUNTSDIR <- "/camp/project/proj-emedlab-vanloo/ICGC_copynumber/battenberg_raw_files/allelecounts/"
# TEMPDIR <- "/camp/project/proj-emedlab-vanloo/home/jdemeul/temp/"
TEMPDIR <- "/tmp/"
ALLELECOUNTSDIR <- LOGRDIR
SNVMNVINDELDIR <- "/camp/project/proj-vanloo-secure/PCAWG/ICGC_snv_mnv/final_consensus_12oct_passonly/"

# BAFDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_rerun_005_input_to_finalConsCopynum_BAFphased/"
BAFDIR <- LOGRDIR

# CNDIR <- "/camp/project/proj-emedlab-vanloo/ICGC_consensus_copynumber/20170119_release/"
# BASEOUT <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/"
BASEOUT <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210702_vafpipeline_out_alphapt1_hetonly/"
# BAMDIR <- "/camp/project/proj-emedlab-vanloo/ICGC/"
PHASINGDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20181021_hetSNPs_all+phasing_out/"
BAMDIR <- "/camp/project/proj-emedlab-vanloo/ICGC/"


### load one time only
releasetable <- read_pcawg_release_table(release_table_file = RELEASETABLEFILE)
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
# load reference alleles
# refalleles_gr <- load_1000G_reference_alleles_new(refallelesdir = REFALLELESDIR)
# saveRDS(object = refalleles_gr, file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/refalleles.RDS", compress = F)
refalleles_gr <- readRDS("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016-17_ICGC/infinite_sites/data/refalleles.RDS")
# immune loci
immune_loci <- GRanges(seqnames = c(14, 7, 7, 14, 2, 22, 6), ranges = IRanges(start = c(22090057, 141998851, 38279625, 106032614, 89156674, 22380474, 28477797),
                                                                              end = c(23021075, 142510972, 38407656, 107288051, 90274235, 23265085, 33448354)), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))

sampleid <- TUMOURNAME

####################
# start with sample ID (i.e. tumour wgs aliquot ID)
# sampleid <- releasetable[1, "tumor_wgs_aliquot_id"]
# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "05780d48-80e7-4d70-b00c-081f8a9519f2"
# sampleid <- "909f3c5d-89fc-419b-a654-75ac1dbb149f"
# sampleid <- "ca8fa9f5-3190-440d-9879-22e33d05ca6c"
# sampleid <- "887616c5-06a7-4e83-948c-3546202349fb"
####################

# run_baflogr_pipeline <- function(sampleid) {
  
  print(paste0("Loading data for sample ", sampleid))
  
  options(bitmapType = "cairo")
  
  sampledir <- file.path(BASEOUT, sampleid)
  dir.create(sampledir)
  #######
  samplerhopsi <- rhopsi[rhopsi$samplename == sampleid, c("purity", "ploidy"), drop = T]
  inferred_sex <- sumtable_whole[sumtable_whole$samplename == sampleid, "inferred_sex", drop = T]
  
  print(paste0("Loading allelecounts for sample ", sampleid))
  allelecounts_tum_gr <- get_tum_allecounts_gr(allelecountsdir = ALLELECOUNTSDIR,
                                               sampleid = sampleid,
                                               bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5,
                                               reference_alleles = refalleles_gr,
                                               tempdir = TEMPDIR)

  if (is.null(allelecounts_tum_gr))# | !dir.exists(sampledir) )
    return(NULL)

  segments_gr <- get_consensus_cn(sampleid = sampleid,
                                  cndir = CNDIR,
                                  bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)

  if (inferred_sex == "male")
    segments_gr <- trim_XY_to_PAR(segments_gr, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)

  print(paste0("Getting phased BAF for sample ", sampleid))
  phasedbaf_gr <- get_phased_BAF(bafdir = BAFDIR,
                                 sampleid = sampleid,
                                 bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5,
                                 allelecounts = allelecounts_tum_gr)

  rm(allelecounts_tum_gr)

  print(paste0("Getting corrected Log R for sample ", sampleid))
  gccorrlogr_gr <- get_GCcorr_logr(logrdir = LOGRDIR,
                                   sampleid = sampleid,
                                   # segments_gr = segments_gr,
                                   bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)

  print(paste0("Loading SNV/MNVs for sample ", sampleid))
  snvs_vcf <- get_snv_mnvs(sampleid = sampleid,
                           snvdir = SNVMNVINDELDIR)
  
  ###### check for load
  if (any(is.null(segments_gr), is.null(phasedbaf_gr), is.null(gccorrlogr_gr), is.null(snvs_vcf))) {
    print(paste0("Missing:", c("CN segments", "Phased BAF", "LogR", "SNVs")[c(is.null(segments_gr), is.null(phasedbaf_gr), is.null(gccorrlogr_gr), is.null(snvs_vcf))], collapse = " "))
    return(NULL)
  }
  ######
  
  ### annotate variants further for QC
  snvs_vcf <- add_snv_qc(sampleid = sampleid,
                               sampledir = sampledir,
                               releasetable = releasetable,
                               snvs_vcf = snvs_vcf,
                               hg19tohg38 = hg19tohg38chain,
                               ncores = NCORES,
                         checkbam = F)
  
  ### additional QC

  print(paste0("Summarising segments for sample ", sampleid))
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
                                  immune_loci = immune_loci,
                                  germline_svs = dbvarcommon1kg)
  
  
  print(paste0("Annotating with phasing and deriving violations for sample ", sampleid))
  # debug(get_metrics)
  # debug(get_performance_metrics)
  # debug(call_parallel_violations)
  finhits <- call_parallel_violations(sampleid = sampleid,
                                      sampledir = sampledir,
                                      phasingdir = PHASINGDIR,
                                      nboot = 100,
                                      alpha = .1)
  
  
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
  
  # return(NULL)
# }


# 
# 
# #### rslurm submission command
# rslurmdf <- releasetable[, "tumor_wgs_aliquot_id", drop = F]
# colnames(rslurmdf) <- "sampleid"

# rslurmdf <- data.frame(sampleid = c("00aa769d-622c-433e-8a8a-63fb5c41ea42", "0bfd1043-5142-3662-e050-11ac0c486501", "0bfd1043-7343-fdd0-e050-11ac0c484cab", 
#                                     "0bfd1043-817c-e3e4-e050-11ac0c4860c5", "0bfd1068-3fc5-a95b-e050-11ac0c4860c3", "0bfd1068-3fe4-a95b-e050-11ac0c4860c3", 
#                                     "0bfe2ac9-0afd-c248-e050-11ac0d487e1c", "10ad692b-4c3d-42de-9b5e-4968441388b3", "154f80bd-984c-4792-bb89-20c4da0c08e0",
#                                     "228fb827-c05e-494c-8a21-e1d925e100cb", "2790b964-63e3-49aa-bf8c-9a00d3448c25", "2df02f2b-9f1c-4249-b3b4-b03079cd97d9",
#                                     "3fba4880-cb7b-4ac5-ab5f-728614faa1ea", "42d20028-0ddc-4dac-9f05-d674f8915f21", "561fd34c-7c7d-4df0-bbfc-3d31147ca562",
#                                     "60c33e32-7e19-4e71-b075-a63fcf27e660", "76a0d9c9-5e69-44e8-9ed2-6d2e387803fc", "8a929c55-35a6-4645-bb70-4b85d281b139",
#                                     "9c27fedd-b1b3-4af0-9e9b-20271854db08", "c95a2b1b-726c-4608-9fff-d57b6f1aa75a", "dc505248-ed04-4f77-a7c6-3fefbc5df27b",
#                                     "e6eda5db-4d4f-418e-b0d4-ed9b3e5259d3", "ea7d37ca-0dac-4ae6-ad03-97c6df3d116d", "f283ed80-8302-4f26-99ed-ea20d101289d",
#                                     "f38b5d2e-5cab-45c7-bb0a-38b2efc5c156"), stringsAsFactors = T)

# existing_outfiles <- list.files(path = BASEOUT, pattern = "_InfSites_VAFpipeline_summarystats\\.txt$", full.names = T, recursive = T)
# # subset to non-empty files
# existing_outfiles <- gsub(pattern = "_InfSites_VAFpipeline_summarystats.txt", replacement = "", 
#                           fixed = T, x = basename(existing_outfiles[which(file.size(existing_outfiles) > 0)]))
# rslurmdf <- rslurmdf[!rslurmdf$sampleid %in% existing_outfiles, , drop = F]
# 
# newrunfiles <- gsub(pattern = "_beagle5", replacement = "", x = list.files(path = LOGRDIR), fixed = T)
# rslurmdf <- rslurmdf[rslurmdf$sampleid %in% newrunfiles, , drop = F]

# debug(run_baflogr_pipeline)
# lapply(X = rslurmdf$sampleid, FUN = run_baflogr_pipeline)

# debug(test_clean_sites)
# debug(annotate_gc_hits)
# debug(call_parallel_violations)
# run_baflogr_pipeline(sampleid = rslurmdf[4,])
# debug(run_baflogr_pipeline)
# run_baflogr_pipeline(sampleid = sampleid)
# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "498ecf81-921d-4df9-a6a4-a625f484e823"
# debug(test_clean_sites)
# debug(run_baflogr_pipeline)
# run_baflogr_pipeline(sampleid = rslurmdf[1,])
# debug(run_baflogr_pipeline)
# for (sname in sample(rslurmdf[,"sampleid"], size = 100, replace = F)) {
#   print(sname)
#   run_baflogr_pipeline(sampleid = sname)
# }
# baflogrjob <- slurm_apply(f = run_baflogr_pipeline, params = rslurmdf[,, drop = F], jobname = "baflogr_run_2021_remainder_3", nodes = min(463, nrow(rslurmdf)), cpus_per_node = 1, add_objects = ls(),
#               pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(exclude = c("worker-himem001")), submit = T)


# baflogrjob <- slurm_apply(f = run_baflogr_pipeline, params = rslurmdf[,, drop = F], jobname = "precrec", nodes = 250, cpus_per_node = 16, add_objects = ls(),
#                           pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(exclude = "fat-worker00[1-4]"), submit = T)


# baflogrjob <- slurm_apply(f = run_baflogr_pipeline, params = rslurmdf[,, drop = F], jobname = "final29", nodes = 7, cpus_per_node = 1, add_objects = ls(),
#                           pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(exclude = "fat-worker00[1-4]"), submit = T)

