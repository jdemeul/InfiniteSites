#### moving and renaming artefact stats files

INDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/code/isa_m2_calls/"

metricfiles <- list.files(path = INDIR, pattern = "_metrics", full.names = T)
donesamples <- gsub(pattern = ".bait_bias_detail_metrics", replacement = "", x = list.files(path = INDIR, pattern = ".bait_bias_detail_metrics", full.names = F), fixed = T)
alldirs <- list.dirs(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/", recursive = T, full.names = T)

for (fle in metricfiles) {
  file.rename(from = metricfiles, to = gsub(pattern = "_.", replacement = ".", x = metricfiles, fixed = T))
}


for (sampleid in donesamples) {
  sampledir <- grep(value = T, pattern = sampleid, x = alldirs)
  system(command = paste0("mv ", file.path(INDIR, sampleid), "* ", sampledir, "/."), wait = T)
}
                    