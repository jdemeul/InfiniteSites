library(parallel)

# aggregate simulator results


# RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

NCORES <- 12


# CNDIR <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"
PERMUTBASE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_all_writefrac/"
SAMPLEBASE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20200603_SampleBasedSim/"

### load one time only
sumtable_whole <- read.delim(file = SUMTABLE_WHOLE, as.is = T)

sampleid <- sumtable_whole[1, "samplename"]
placeholder <- read.delim(file = file.path(PERMUTBASE, "0009b464-b376-4fbc-8a56-da538269a02f", "0009b464-b376-4fbc-8a56-da538269a02f_infsites_permut_effgenfrac1.txt"), as.is = T)[, c(1:3)]
# placeholder$counts_tot <- NA

read_sim_results <- function(sampleid, placeholder) {
  permbasedresfile <- file.path(PERMUTBASE, sampleid, paste0(sampleid, "_infsites_permut_effgenfrac1.txt"))
  samplebasedresfile <- file.path(SAMPLEBASE, sampleid, paste0(sampleid, "_infsites_permut_samplebasedsim.txt"))
  
  
  if (file.exists(permbasedresfile)) {
    permres <- read.delim(file = permbasedresfile, as.is = T)[,7]
  } else {
    permres <- rep(NA, 456)
  }
  
  if (file.exists(samplebasedresfile)) {
    sampleres <- read.delim(file = samplebasedresfile, as.is = T)[,7]
  } else {
    sampleres <- rep(NA, 456)
  }
  
  return(list(perm = permres, sample = sampleres))
}




allsimdata <- mclapply(X = sumtable_whole$samplename, FUN = read_sim_results, placeholder = placeholder, mc.preschedule = T, mc.cores = 12)
permutdata <- cbind(placeholder, do.call(cbind, lapply(X = allsimdata, FUN = function(x) x[["perm"]])))
colnames(permutdata)[-c(1:3)] <- sumtable_whole$samplename
sampledata <- cbind(placeholder, do.call(cbind, lapply(X = allsimdata, FUN = function(x) x[["sample"]])))
colnames(sampledata)[-c(1:3)] <- sumtable_whole$samplename

write.table(x = permutdata, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210304_permutation_based_sim_totals1000x_compiled.txt", quote = F, sep = "\t", row.names = F)
write.table(x = sampledata, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210304_sample_based_sim_totals1000x_compiled.txt", quote = F, sep = "\t", row.names = F)
