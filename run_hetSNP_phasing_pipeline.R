### Run heterozygous SNP collection and phasing pipelines

library(readr)
library(GenomicRanges)
library(VariantAnnotation)
library(Rsamtools)
library(parallel)
library(rslurm)

# genome
library(BSgenome.Hsapiens.1000genomes.hs37d5)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/getCleanHetSNPs.R")
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/phasing.R")
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")

# SAMPLE <- "93ff786e-0165-4b02-8d27-806d422e93fc"
# SAMPLE <- "deb9fbb6-656b-41ce-8299-554efc2379bd"

max_distance <- 700
minMapQ <- 20
minBaseQ <- 25
max_mm_read <- 2
max_mm_pair <- 3
NCORES <- 10

REFALLELESDIR <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_loci/"
PROBLOCIFILE <- "/srv/shared/vanloo/pipeline-files/human/references/battenberg/battenberg_problem_loci/probloci_270415.txt.gz"
RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
ALLELECOUNTSTEMPDIR <- "/home/jdemeul/temp/"

# ALLELECOUNTSDIR <- "/srv/shared/vanloo/ICGC_copynumber/battenberg_raw_files/combined_allelecounts/"
ALLELECOUNTSDIR <- "/srv/shared/vanloo/ICGC_copynumber/battenberg_raw_files/allelecounts/"
SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/"
CNDIR <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"
BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20181021_hetSNPs_all+phasing_out/"
BAMDIR <- "/srv/shared/vanloo/ICGC/"


### load one time only
releasetable <- read_pcawg_release_table(release_table_file = RELEASETABLEFILE)
genome_seqinfo <- GenomeInfoDb::keepSeqlevels(x = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5), value = c(1:22,"X","Y"))
refalleles <- load_reference_alleles(refallelesdir = REFALLELESDIR, chrominfo = genome_seqinfo)
probloci <- load_problematic_loci(problocifile = PROBLOCIFILE, chrominfo = genome_seqinfo)


# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"

### get paths to all required files
# tumour_bam_file <- "/srv/shared/vanloo/ICGC/MELA-AU/WGS/84e9d222d5b28a90b27643850e6eb513.bam"
# normal_bam_file <- "/srv/shared/vanloo/ICGC/MELA-AU/WGS/0a69ff126ac2404f505c75a18cb7f131.bam"

run_phasing_pipeline <- function(sampleid) {

  print(paste0("Running sample ", sampleid))
  outdir <- file.path(BASEOUT, sampleid)
  dir.create(outdir, showWarnings = F)
  
  # if (length(list.files(outdir)) == 3) {
  #   return(NULL)
  # }
  
  indelfile <- list.files(path = SNVMNVINDELDIR, pattern = paste0(sampleid, ".consensus.20161006.somatic.indel.vcf.gz$"), full.names = T, recursive = T)
  breaksfile <- list.files(path = CNDIR, pattern = paste0(sampleid, ".consensus.20170119.somatic.cna.annotated.txt$"), full.names = T)
  snv_mnvfile <- list.files(path = SNVMNVINDELDIR, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz$"), full.names = T, recursive = T)
  # allelecountsfile <- file.path(ALLELECOUNTSDIR, paste0(sampleid, "_alleleCounts.tab.gz"))
  allelecountsfile_normal <- file.path(ALLELECOUNTSDIR, paste0(releasetable[releasetable$tumor_wgs_aliquot_id == sampleid, "normal_wgs_aliquot_id"], "_allelecounts.tar.gz"))
  
  germ_loci_file <- file.path(outdir, paste0(sampleid, "_germ_loci.txt.gz"))
  
  if (file.exists(allelecountsfile_normal) & !file.exists(germ_loci_file)) {
    print(paste0("Getting hetSNPs"))
    get_clean_hetSNPs(sampleid = sampleid,
                      allelecountsfile = allelecountsfile_normal,
                      indelfile = indelfile,
                      breaksfile = breaksfile,
                      probloci = probloci,
                      chrominfo = genome_seqinfo,
                      HDI90 = c(lower = 1/3, upper = 2/3),
                      mindepth = 10, # modified
                      heterozygousfilter = .1, # modified
                      alleles = refalleles,
                      outfile = file.path(outdir, paste0(sampleid, "_germ_loci.txt")),
                      indelfilter = F, # last three added in rerun to include all (potentially) phaseable loci
                      breakfilter = F,
                      hdifilter = F)
    system(command = paste0("gzip ", germ_loci_file), wait = T)
  }

  
  tumour_bam_file <- list.files(path = BAMDIR, pattern = paste0(releasetable[releasetable$tumor_wgs_aliquot_id == sampleid, "tumor_wgs_bwa_alignment_bam_file_name"], "$"), full.names = T, recursive = T)
  # normal_bam_file <- list.files(path = BAMDIR, pattern = paste0(releasetable[releasetable$tumor_wgs_aliquot_id == sampleid, "normal_wgs_bwa_alignment_bam_file_name"], "$"), full.names = T, recursive = T)
  
  if (length(tumour_bam_file) == 0 | all(!file.exists(allelecountsfile_normal, snv_mnvfile))) {
    return(NULL)
  }

  print(paste0("Phasing tumour"))
  snv_phasing(germ_loci_file = germ_loci_file,
              som_vcf_file = snv_mnvfile,
              outfile = file.path(outdir, paste0(sampleid, "_tumour_snv-snp_phased.txt")),
              bam_file = tumour_bam_file,
              bai_file = paste0(tumour_bam_file, ".bai"),
              max_distance = max_distance,
              minMapQ = minMapQ,
              max_mm_read = max_mm_read,
              max_mm_pair = max_mm_pair,
              minBaseQ = minBaseQ,
              ncores = NCORES, preschedule = F,
              chrominfo = genome_seqinfo)
  
  system(command = paste0("gzip -f ", file.path(outdir, paste0(sampleid, "_tumour_snv-snp_phased.txt"))), wait = T)
  
  if (length(normal_bam_file) == 0) {
    return(NULL)
  }
  
  print(paste0("Phasing normal"))
  snv_phasing(germ_loci_file = file.path(outdir, paste0(sampleid, "_germ_loci.txt")),
              som_vcf_file = snv_mnvfile,
              outfile = file.path(outdir, paste0(sampleid, "_normal_snv-snp_phased.txt")),
              bam_file = normal_bam_file,
              bai_file = paste0(normal_bam_file, ".bai"),
              max_distance = max_distance,
              minMapQ = minMapQ,
              max_mm_read = max_mm_read,
              max_mm_pair = max_mm_pair,
              minBaseQ = minBaseQ,
              ncores = NCORES, preschedule = T,
              chrominfo = genome_seqinfo)
  
  return(NULL)
}

# passedfiles <- names(which(table(dirname(list.files(path = BASEOUT, recursive = T, all.files = F))) == 3))
# rslurmdf <- releasetable[!releasetable$tumor_wgs_aliquot_id %in% passedfiles, "tumor_wgs_aliquot_id" , drop = F]
# missingbams <- list.files(path = "/srv/shared/vanloo/ICGC_missingbam/", pattern = "*.bam$")
# rslurmdf <- releasetable[releasetable$tumor_wgs_bwa_alignment_bam_file_name %in% missingbams, "tumor_wgs_aliquot_id" , drop = F]

tumourphasedsamples <- gsub(pattern = "_tumour_snv-snp_phased.txt.gz", replacement = "", x = basename(list.files(path = BASEOUT, pattern = "*_tumour_snv-snp_phased.txt.gz", recursive = T)))
rslurmdf <- releasetable[which(!releasetable$tumor_wgs_aliquot_id %in% tumourphasedsamples), "tumor_wgs_aliquot_id" , drop = F]

# rslurmdf <- releasetable[, "tumor_wgs_aliquot_id" , drop = F]
colnames(rslurmdf) <- "sampleid"

# normalphasedsamples <- gsub(pattern = "_normal_snv-snp_phased.txt", replacement = "", x = basename(list.files(path = BASEOUT, pattern = "*_normal_snv-snp_phased.txt$", full.names = F, recursive = T)))
# tumourphasedsamples <- gsub(pattern = "_tumour_snv-snp_phased.txt", replacement = "", x = basename(list.files(path = BASEOUT, pattern = "*_tumour_snv-snp_phased.txt$", full.names = F, recursive = T)))

#### some files may not be complete ... check whether all chroms present
# tumourphasedsamples <- list.files(path = BASEOUT, pattern = "*_tumour_snv-snp_phased.txt.gz", full.names = T, recursive = T)
# completesamples <- mclapply(X = tumourphasedsamples, FUN = function(x) all(c(1:22, "X") %in% unique(read_tsv(file = x, col_types = c("c_____________")))$chr), mc.preschedule = T, mc.cores = NCORES)
# samplestorerun <- gsub(pattern = "_tumour_snv-snp_phased.txt", replacement = "", x = basename(tumourphasedsamples))[!unlist(completesamples)]
# lapply(X = samplestorerun, FUN = run_phasing_pipeline)



# rslurmdf <- rslurmdf[which(!(rslurmdf$sampleid %in% tumourphasedsamples)), ]
# lapply(X = rslurmdf, FUN = run_phasing_pipeline)

# SAMPLE <- "6c8f3dc9-21bf-4859-9599-231ac040eb7d"
# debug(get_clean_hetSNPs)
# debug(snv_phasing)
# debug(snv_phasing_chr)
# debug(load_allelecounts_normal)
# debug(run_phasing_pipeline)
# run_phasing_pipeline(sampleid = SAMPLE)
# sampleidtest <- rslurmdf[1,]

phasingjob <- slurm_apply(f = run_phasing_pipeline, params = rslurmdf[,, drop = F], jobname = "phase_missing_run_all", nodes = 25, cpus_per_node = 1, add_objects = ls(),
            pkgs = rev(.packages()), libPaths = .libPaths(), submit = T, slurm_options = list(exclude = "fat-worker00[1-6]"))

# print_job_status(phasingjob)
# cancel_slurm(phasingjob)

