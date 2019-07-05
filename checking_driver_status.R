# #### check PCAWG variant calls again vs the double parallel ones
# 
# vafhitsfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/20180529_hetSNPs+phasing_out/", pattern = "_baflogr_gc_results.txt", full.names = T, recursive = T)
# vafhits <- mclapply(X = vafhitsfiles, FUN = read.delim, as.is = T, mc.cores = 16)
# vafhits <- lapply(X = vafhits, FUN = function(x) x[x$type == "snv" & x$padj < .05 & x$hetsnpneighbours > 0 & x$nsnphitsoutofrange/x$hetsnpneighbours < .1 & x$nlogrhits/x$nlogrneighbours < .1,  ])
# 
# vafhitsdf <- do.call(rbind, vafhits)
# vafhitsdf$sampleid <- rep(x = gsub(pattern = "_baflogr_gc_results.txt", replacement = "", x = basename(vafhitsfiles)), sapply(X = vafhits, FUN = nrow))
# 
# driversdf <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
# 
# parallelhits <- vafhitsdf[paste0(vafhitsdf$sampleid, "_", vafhitsdf$seqnames, "_", vafhitsdf$start) %in% paste0(driversdf$sample, "_", driversdf$chr, "_", driversdf$pos), ]
# 


library(metap)
### parallel clean hit samples
parhitsamples <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2.txt", as.is = T)

# read actual hits
parhitsfiles <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", parhitsamples$sampleid, "/", parhitsamples$sampleid, "_snv_mnv_infSites_finalhits.txt")
parhits <- lapply(X = parhitsfiles, FUN = read.delim, as.is = T)
parhitsdf <- do.call(rbind, parhits)
parhitsdf$sampleid <- rep(parhitsamples$sampleid, sapply(parhits, nrow))

parhitsdf$bafpval_comb <- rep(1, nrow(parhitsdf))
goodbafidxs <- which(!(is.na(parhitsdf$bafpval_pre) | is.na(parhitsdf$bafpval_post)))
parhitsdf$bafpval_comb[goodbafidxs] <- apply(X = parhitsdf[goodbafidxs, c("bafpval_pre", "bafpval_post")], MARGIN = 1, FUN = function(x) sumlog(p = x)$p)

parhitsdf$logrpval_comb <- rep(1, nrow(parhitsdf))
goodlogridxs <- which(!(is.na(parhitsdf$logrpval_pre) | is.na(parhitsdf$logrpval_post)))
parhitsdf$logrpval_comb[goodlogridxs] <- apply(X = parhitsdf[goodlogridxs, c("logrpval_pre", "logrpval_post")], MARGIN = 1, FUN = function(x) sumlog(p = x)$p)

parhitsdf <- parhitsdf[parhitsdf$bafpval_comb > 1e-2 & parhitsdf$logrpval_comb > 1e-2, ]

parhitsdf$ref <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$ref))
parhitsdf$alt <- gsub(pattern = "RUE", replacement = "", x = as.character(parhitsdf$alt))

write.table(x = parhitsdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt", quote = F, row.names = F, col.names = T, sep = "\t")

parhitsdf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt", as.is = T, sep = "\t")

driversdf <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
# pardrivers <- parhitsdf[paste0(parhitsdf$sampleid, "_", parhitsdf$start) %in% paste0(driversdf$sample, "_", driversdf$pos), ]
pardrivers <- parhitsdf[paste0(parhitsdf$sampleid, "_", parhitsdf$start) %in% paste0(driversdf$sample, "_", driversdf$pos), ]




#### check reannotated drivers
pardriv <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/annovar/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.hg19_multianno.txt", as.is = T)
pardrivorig <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt", as.is = T, header = T)

all(pardriv$Start == pardrivorig$start)

cgc <- read.csv(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/cancer_gene_census.csv", as.is = T)
driversdf <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS1_compendium_mutational_drivers_20180110.txt", as.is = T)

sumtab <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt", as.is = T)

cgidxs <- which(sapply(X = pardriv$Gene.refGene, FUN = function(x, y) any(strsplit(x = x, split = ";", fixed = T)[[1]] %in% y), y = unique(cgc$Gene.Symbol, driversdf$Gene)))
# cgidxs <- 1:nrow(pardrivorig)
pardrivhit <- cbind(pardrivorig[cgidxs, ], pardriv[cgidxs, -c(1:5)])

pardrivhit$histology_abbreviation <- sumtab[match(x = pardrivhit$sampleid, table = sumtab$samplename), "histology_abbreviation"]



thirddriv1 <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/annovar/InfSitesBiallelicM2recall_variants_alt1.hg19_multianno.txt", as.is = T)
thirddriv2 <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/annovar/InfSitesBiallelicM2recall_variants_alt2.hg19_multianno.txt", as.is = T)
thirddrivorig <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T)

cgidxsthird <- which(sapply(X = thirddriv1$Gene.refGene, FUN = function(x, y) any(strsplit(x = x, split = ";", fixed = T)[[1]] %in% y), y = unique(cgc$Gene.Symbol, driversdf$Gene)) |
                        sapply(X = thirddriv2$Gene.refGene, FUN = function(x, y) any(strsplit(x = x, split = ";", fixed = T)[[1]] %in% y), y = unique(cgc$Gene.Symbol, driversdf$Gene)))

thirddrivhit <- cbind(thirddrivorig[cgidxsthird, ], thirddriv1[cgidxsthird, -c(1:5)], thirddriv2[cgidxsthird, -c(1:5)])

thirddrivhit$histology_abbreviation <- sumtab[match(x = thirddrivhit$sampleid, table = sumtab$samplename), "histology_abbreviation"]

pardrivhitex <- pardrivhit[pardrivhit$Func.refGene == "exonic", ]
thirddrivhitex <- thirddrivhit[thirddrivhit$Func.refGene == "exonic", ]

write.table(x = pardrivhitex, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/annovar/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_exondrivers.txt.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = thirddrivhitex, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/annovar/InfSitesBiallelicM2recall_variants_v2_exondrivers.txt", quote = F, sep = "\t", row.names = F, col.names = T)



