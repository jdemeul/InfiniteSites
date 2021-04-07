### annotate vcf with normal bam file depth

library(rslurm)

SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/"
RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
PARHITSRES <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt"
DIVHITSRES <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary.txt"
BAMDIR <- "/srv/shared/vanloo/ICGC/"
OUTDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/normaldepthannotatedvcfs/"

sumtab <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
releasetab <- read.delim(file = RELEASETABLEFILE, as.is = T)
samplesofinterest <- unique(c(read.delim(file = PARHITSRES, as.is = T)$sampleid, read.delim(file = DIVHITSRES, as.is = T)$sampleid))



annotate_vcf_depth <- function(sampleid) {
  print(paste0("Running sample ", sampleid))
  
  snv_mnvfile <- list.files(path = SNVMNVINDELDIR, recursive = T, full.names = T, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz$"))
  normalbamfile <- list.files(path = BAMDIR, recursive = T, full.names = T, pattern = paste0(releasetab$normal_wgs_bwa_alignment_bam_file_name, "$"))
  outfile <- file.path(OUTDIR, paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.wdepth.vcf"))
  
  if (all(file.exists(snv_mnvfile, normalbamfile)) & !file.exists(outfile)) {
    systemcmd <- paste0("/srv/sw/eb/software/GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8/gatk  --java-options '-Xmx4G' AnnotateVcfWithBamDepth",
                        " -V ", snv_mnvfile,
                        " -I ", normalbamfile,
                        " -O ", outfile)
    print(x = systemcmd)
    system(command = systemcmd)
  }
  
  return(NULL)
}


# rslurmdf <- data.frame(sampleid = samplesofinterest)
rslurmdf <- data.frame(sampleid = sumtab$samplename)
# annotate_vcf_depth(samplesofinterest[1])

addnormaldepth <- slurm_apply(f = annotate_vcf_depth, params = rslurmdf[,, drop = F], jobname = "adddepth", nodes = length(samplesofinterest)/5, cpus_per_node = 5, add_objects = ls(),
                              pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(exclude = "fat-worker00[1-4]"), submit = T)
