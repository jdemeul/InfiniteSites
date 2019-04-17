### check phaseing of VAF hit SNVs (moving to InfSites pipeline rather than GC)

library(ggplot2)
library(ROCR)
# library(boot)

VAFHITSDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out/"
PHASINGDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20181021_hetSNPs_all+phasing_out/"


### fions

get_performance_metrics <- function(df, plotting = F, sampleid = sampleid, sampledir = OUTDIR) {
  
  if (sum(df$is_confirmed) == 0) {
    return(c(sensitivity = NA, specificity = NA, cutoff = NA, auc = NA))
  }
  
  infsitespred <- prediction(predictions = -log10(df$pval), labels = df$is_confirmed)
  infsitesperf <- performance(prediction.obj = infsitespred, measure = "tpr", x.measure = "fpr")
  
  optimperf <- opt_cut(perf = infsitesperf, pred = infsitespred)
  optauc <- performance(prediction.obj = infsitespred, measure = "auc")@y.values[[1]]
  
  if (plotting) {
    p6 <- ggplot(data = as.data.frame(t(optimperf)), mapping = aes(y = sensitivity, x = 1 - specificity))
    p6 <- p6 + geom_label(mapping = aes(label = paste0("-log10(pval) = ", round(cutoff, digits = 2))), nudge_x = .05, nudge_y = -.05, hjust = "left", vjust = "top") + geom_point(color = "red")
    p6 <- p6 + geom_line(data = data.frame(TPR = infsitesperf@y.values[[1]], FPR = infsitesperf@x.values[[1]]), mapping = aes(x = FPR, y = TPR))
    p6 <- p6 + theme_minimal() + coord_equal(xlim = c(0,1), ylim = c(0,1)) + geom_abline(color = "grey", linetype = "dashed")
    # p6 <- p6 + labs(title = paste0("sample ", sampleid, " - ", "\ntotal testable: ", nrow(vafhitsdf), " - parallel violations: ", nrow(finalhits), " / total phaseable: ", nrow(vafhitsroc), " - violation confirmed: ", sum(vafhitsroc$is_confirmed)))
    
    ggsave(filename = file.path(sampledir, paste0(sampleid, "_ROC.png")), plot = p6, width = 10, height = 10, units = "cm")
  }
  return(c(optimperf["sensitivity", 1], optimperf["specificity", 1], 10^-optimperf["cutoff", 1], auc = optauc))
}


opt_cut <- function(perf, pred){
  cut_ind = mapply(FUN=function(x, y, p) {
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, x = perf@x.values, y = perf@y.values, p = pred@cutoffs)
  return(cut_ind)
}
#### fions


# sampleid <- "66eb4833-1b87-4fd9-a53d-26dc7ad6de29"
sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "93ff786e-0165-4b02-8d27-806d422e93fc"

sampledir <- file.path(VAFHITSDIR, sampleid)

# read data if exists, otherwise exit
vafhitsfile <- file.path(VAFHITSDIR, sampleid, paste0(sampleid, "_snv_mnv_infSites_annotated.txt"))
phasinghitsfile <- file.path(PHASINGDIR, sampleid, paste0(sampleid, "_tumour_snv-snp_phased.txt"))

if (any(!file.exists(vafhitsfile, phasinghitsfile))) return(NULL)

vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
phasinghits <- read.delim(file = phasinghitsfile, as.is = T)

# SNV-SNP phasing only directly informative
phasinghits <- phasinghits[which(phasinghits$type1 != phasinghits$type2), ]

# make sure there is a decent amount of pair coverage necessary to see anything (i.e. 2 reads for each SNP allele, 4 reads from the somatic variant allele, fair since we won't pick up anything subclonal anyhow)
coveredphasinghits <- which(ifelse(phasinghits$type1 == "SNP",
                             phasinghits$Num_ref_ref + phasinghits$Num_ref_alt >= 2 & phasinghits$Num_alt_ref + phasinghits$Num_alt_alt >= 2 & phasinghits$Num_ref_alt + phasinghits$Num_alt_alt >= 4 ,
                             phasinghits$Num_ref_ref + phasinghits$Num_alt_ref >= 2 & phasinghits$Num_alt_alt + phasinghits$Num_ref_alt >= 2 & phasinghits$Num_alt_ref + phasinghits$Num_alt_alt >= 4))
phasinghits <- phasinghits[coveredphasinghits, ]

# create snv indices for quick overlaps
phasinghits$snvidx <- ifelse(phasinghits$type1 == "SNV", paste0(phasinghits$chr, "_", phasinghits$pos1),  paste0(phasinghits$chr, "_", phasinghits$pos2))
vafhitsidxs <- paste0(vafhitsdf$chr, "_", vafhitsdf$start)

# subset to the ones which have been evaluated for ISA violations
phasinghitssub <- phasinghits[which(phasinghits$snvidx %in% vafhitsidxs), ]
# and call violations from phasing data
phasinghitssub$totalpwcov <- rowSums(phasinghitssub[, c("Num_ref_ref", "Num_ref_alt", "Num_alt_alt", "Num_alt_ref")])
phasinghitssub$isaviol <- ifelse(phasinghitssub$type1 == "SNP",
                                 phasinghitssub$Num_ref_alt > 1 & phasinghitssub$Num_alt_alt > 1 & phasinghitssub$Num_ref_alt/phasinghitssub$totalpwcov > .1 & phasinghitssub$Num_alt_alt/phasinghitssub$totalpwcov > .1,
                                 phasinghitssub$Num_alt_ref > 1 & phasinghitssub$Num_alt_alt > 1 & phasinghitssub$Num_alt_ref/phasinghitssub$totalpwcov > .1 & phasinghitssub$Num_alt_alt/phasinghitssub$totalpwcov > .1)

# assign calls
vafhitsdf$is_phaseable <- vafhitsidxs %in% phasinghitssub$snvidx
vafhitsdf$is_confirmed <- vafhitsidxs %in% names(which(c(by(data = phasinghitssub$isaviol, INDICES = phasinghitssub$snvidx, FUN = any))))

# write out the phasing-annotated variants and associated phasingdata

# prepare threshold optimisation
# vafhitsdfsub <- vafhitsdf[, ]

# plotting during devel
# p1 <- ggplot(data = vafhitsdfsub, mapping = aes(y = pval, x = is_confirmed)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + scale_y_log10()
# p1
# p2 <- ggplot(data = vafhitsdfsub, mapping = aes(y = pval, x = is_confirmed)) + geom_jitter(mapping = aes(colour = nhetsnps25bp >= 2)) + scale_y_log10()
# p2
# p2 <- ggplot(data = vafhitsdfsub[which(vafhitsdfsub$pval < 1e-4), ], mapping = aes(y = pmin(bafpval_pre, bafpval_post), x = is_confirmed)) + geom_jitter(mapping = aes(colour = pmin(bafpval_pre, bafpval_post) > 1e-3)) + scale_y_log10()
# p2
# p3 <- ggplot(data = vafhitsdfsub[which(vafhitsdfsub$pval < 1e-4), ], mapping = aes(y = pmin(logrpval_pre, logrpval_post), x = is_confirmed)) + geom_jitter(mapping = aes(colour = pmin(logrpval_pre, logrpval_post) > 1e-3)) + scale_y_log10()
# p3
# p3 <- ggplot(data = vafhitsdfsub[which(vafhitsdfsub$pval < 1e-4), ], mapping = aes(y = pval, x = is_confirmed, colour = is.na(total_cn) | minor_cn == 0)) + geom_jitter() + scale_y_log10()
# p3
# p4 <- ggplot(data = vafhitsdfsub, mapping = aes(y = pval, x = is_confirmed)) + geom_boxplot(outlier.shape = NA) + geom_jitter(mapping = aes(colour = pmin(bafpval_pre, bafpval_post) > 1e-3 & pmin(logrpval_pre, logrpval_post) > 1e-3)) + scale_y_log10()
# p4

# filter out variants where adjacent (het)SNPs do not behave (BAF/LogR) according to segment
# filter out variants in the immune regions (as allele frequencies may be messed up)
# also filter out sites which have ≥ 2 hetSNPs within 25bp window as they considerably bias the allelecounts when phased
vafhitsroc <- vafhitsdf[which(vafhitsdf$is_phaseable & !is.na(vafhitsdf$pval) &
                                pmin(vafhitsdf$bafpval_pre, vafhitsdf$bafpval_post) > 1e-3 &
                                   pmin(vafhitsdf$logrpval_pre, vafhitsdf$logrpval_post) > 1e-3 &
                                   !vafhitsdf$immune_locus & vafhitsdf$nhetsnps25bp < 2),  ]

#summary stats
par_phasing_conf <- sum(vafhitsroc$is_confirmed)

if (par_phasing_conf > 0) {
  bootout <- boot(data = vafhitsroc, statistic = function(x, i) {get_performance_metrics(df = x[i,])}, R = 1000)
  perfmets <- setNames(object = colMeans(bootout$t, na.rm = T), names(bootout$t0))
} else if (nrow(vafhitsdf) > 0 ){
  perfmets <- c(sensitivity = NA, specificity = NA, cutoff = sum(p.adjust(vafhitsdf$pval, method = "fdr") <= .05 , na.rm = T)*.05/sum(!is.na(vafhitsdf$pval)), auc = NA)
} else {
  perfmets <- c(sensitivity = NA, specificity = NA, cutoff = NA, auc = NA)
}


finalhits <- vafhitsdf[which(vafhitsdf$pval <= perfmets[["cutoff"]] & pmin(vafhitsdf$bafpval_pre, vafhitsdf$bafpval_post) > 1e-3 & pmin(vafhitsdf$logrpval_pre, vafhitsdf$logrpval_post) > 1e-3 & !vafhitsdf$immune_locus & vafhitsdf$nhetsnps25bp < 2), ]
sumstats <- c(tot_testable = nrow(vafhitsdf), nparallel = nrow(finalhits), tot_phaseable = nrow(vafhitsroc), npar_phased = par_phasing_conf, perfmets)

# writing output
write.table(x = finalhits, file = file.path(sampledir, paste0(sampleid, "_snv_mnv_infSites_finalhits.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = vafhitsdf, file = file.path(sampledir, paste0(sampleid, "_snv_mnv_infSites+phasing_annotated.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = phasinghitssub, file = file.path(sampledir, paste0(sampleid, "_tumour_snv-snp_phased_InfSitesInformative.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = sumstats, file = file.path(sampledir, paste0(sampleid, "_InfSites_VAFpipeline_summarystats.txt")), sep = "\t", quote = F, col.names = F, row.names = T)




# 
# #### quick check of SNV slopes in QQ plots
# VAFOUTDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/"
# slopefiles <- list.files(path = VAFOUTDIR, pattern = "_pseudocount_calibrated.txt", full.names = T, recursive = T)
# slopes <- do.call(rbind, lapply(X = slopefiles, FUN = read.delim, as.is = T, header = F))
# slopes$sampleid <- gsub(pattern = "_pseudocount_calibrated.txt", replacement = "", x = basename(slopefiles))
