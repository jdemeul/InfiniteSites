### check phaseing of VAF hit SNVs

sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
sampleid <- "93ff786e-0165-4b02-8d27-806d422e93fc"
library(ggplot2)

vafhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/", sampleid, paste0(sampleid, "_baflogr_gc_results.txt"))
phasinghitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20181021_hetSNPs_all+phasing_out/", sampleid, paste0(sampleid, "_tumour_snv-snp_phased.txt"))

vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
vafhitsdf <- vafhitsdf[vafhitsdf$type == "snv" & vafhitsdf$hetsnpneighbours > 0 & vafhitsdf$nsnphitsoutofrange/vafhitsdf$hetsnpneighbours < .1 & vafhitsdf$nlogrhits/vafhitsdf$nlogrneighbours < .1,  ]
#  & vafhitsdf$padj < .01 
phasinghits <- read.delim(file = phasinghitsfile, as.is = T)
phasinghits <- phasinghits[which(phasinghits$type1 != phasinghits$type2), ]
coveredphasinghits <- which(ifelse(phasinghits$type1 == "SNP",
                             phasinghits$Num_ref_ref + phasinghits$Num_ref_alt > 0 & phasinghits$Num_alt_ref + phasinghits$Num_alt_alt > 0 & phasinghits$Num_ref_alt + phasinghits$Num_alt_alt > 1,
                             phasinghits$Num_ref_ref + phasinghits$Num_alt_ref > 0 & phasinghits$Num_alt_alt + phasinghits$Num_ref_alt > 0 & phasinghits$Num_alt_ref + phasinghits$Num_alt_alt > 1))
phasinghits <- phasinghits[coveredphasinghits, ]

phasidxs1 <- paste0(phasinghits$chr, "_", phasinghits$pos1)
phasidxs2 <- paste0(phasinghits$chr, "_", phasinghits$pos2)

vafhitsidxs <- paste0(vafhitsdf$seqnames, "_", vafhitsdf$start)

phasinghitssub <- phasinghits[which(phasidxs1 %in% vafhitsidxs | phasidxs2 %in% vafhitsidxs), ]
phasinghitssub$isaviol <- ifelse(phasinghitssub$type1 == "SNP",
                                 phasinghitssub$Num_ref_alt > 0 & phasinghitssub$Num_alt_alt > 0,
                                 phasinghitssub$Num_alt_ref > 0 & phasinghitssub$Num_alt_alt > 0)
phasinghitssub$snvidx <- ifelse(phasinghitssub$type1 == "SNP", paste0(phasinghitssub$chr, "_", phasinghitssub$pos2),  paste0(phasinghitssub$chr, "_", phasinghitssub$pos1))

vafhitsdf$is_phaseable <- vafhitsidxs %in% c(phasidxs1, phasidxs2)
vafhitsdf$is_confirmed <- vafhitsidxs %in% names(which(c(by(data = phasinghitssub$isaviol, INDICES = phasinghitssub$snvidx, FUN = any))))
vafhitsdf$logpadj <- -log10(vafhitsdf$padj)

vafhitsdfsub <- vafhitsdf[vafhitsdf$is_phaseable, ]

pairs(vafhitsdf[vafhitsdf$is_phaseable, c("is_confirmed", "logpadj", "hetsnpneighbours", "nsnphitsoutofrange", "nlogrhits")])

sum(vafhitsdf$is_phaseable)
sum(vafhitsdf$is_confirmed)

# library(ROCR)
# mean(sample(x = vafhitsdf[vafhitsdf$is_confirmed, "padj"], size = 1000, replace=T) < sample(x = vafhitsdf[!vafhitsdf$is_confirmed, "padj"], size = 1000, replace=T))

temp <- prediction(predictions = -log10(vafhitsdfsub$padj), labels = vafhitsdfsub$is_confirmed)
temp2 <- performance(prediction.obj = temp, measure = "tpr", x.measure = "fpr")
plot(temp2)

opt.cut <- function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p) {
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, x = perf@x.values, y = perf@y.values, p = pred@cutoffs)
}

print(opt.cut(perf = temp2, pred = temp))
