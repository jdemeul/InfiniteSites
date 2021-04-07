### Rheinbay-type analysis


# load all variants
library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(Biostrings)
library(VariantAnnotation)
library(parallel)
library(ggplot2)

SNVMNVINDELDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/"
SUMTABLE_WHOLE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"

# sampleid <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
# sampleid <- "fc876f5c-8339-bc9c-e040-11ac0d485160"
# sampleid <- "7cae6c0b-36fe-411b-bbba-093a4c846d84"
# sampleid <- "0332b017-17d5-4083-8fc4-9d6f8fdbbbde"
# sampleid <- "804ffa2e-158b-447d-945c-707684134c87"

sumtab <- read.delim(file = SUMTABLE_WHOLE, as.is = T)
sampleids <- sumtab[sumtab$is_preferred, "samplename"]

# parhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2.txt", as.is = T), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
parhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_cleanhits_variants_v2_2021.txt", as.is = T), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
parhits <- parhits[!seqnames(parhits) == "Y"]
parhits <- parhits[parhits$sampleid %in% sampleids]
# thirdhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt", as.is = T), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
thirdhits <- GRanges(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_20210217.txt", as.is = T), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
thirdhits <- thirdhits[!seqnames(thirdhits) == "Y"]
thirdhits <- thirdhits[thirdhits$sampleid %in% sampleids]
thirdhits_ext <- thirdhits[rep(1:length(thirdhits), ifelse(thirdhits$in_PCAWG == "", 2, 1))]
# thirdhits_ext$alt1 <- ifelse(duplicated(paste0(seqnames(thirdhits_ext), ":", start(thirdhits_ext), "_", thirdhits_ext$smapleid)), thirdhits_ext$alt2, thirdhits_ext$alt1)

vcffiles <- paste0(SNVMNVINDELDIR, "/", sampleids, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
vcffilesalt <- paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/graylist/snv_mnv/", sampleids, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
sum(file.exists(vcffiles)) + sum(file.exists(vcffilesalt))
vcffiles <- ifelse(file.exists(vcffiles), vcffiles, vcffilesalt)


allsnvsgr <- sort(unlist(GRangesList(mclapply(X = vcffiles, FUN = function(x) {
  y <- readVcf(file = x, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5),
               param = ScanVcfParam(fixed = c("ALT"), info = c("Variant_Classification")))
  z <- rowRanges(y)
  mcols(z) <- cbind(mcols(z),
                    sampleid = rep(gsub(pattern = ".consensus.20160830.somatic.snv_mnv.vcf.gz", replacement = "", x = basename(x), fixed = T), length(z)),
                    info(y))
  return(z)
  }, mc.preschedule = T, mc.cores = 16))))

allsnvpos <- sort(unique(c(granges(parhits), granges(thirdhits_ext), granges(allsnvsgr))))
allsnvpos$count <- Rle(countOverlaps(query = allsnvpos, subject = allsnvsgr, type = "equal"))
allsnvpos$count_par <- Rle(countOverlaps(query = allsnvpos, subject = parhits, type = "equal"))
allsnvpos$count_div <- Rle(countOverlaps(query = allsnvpos, subject = thirdhits_ext, type = "equal"))
allsnvpos$total <- Rle(allsnvpos$count + allsnvpos$count_par + allsnvpos$count_div)
allsnvpos$bial <- Rle(allsnvpos$count_par+allsnvpos$count_div > 0)
  
# table(allsnvpos$count+allsnvpos$count_par+allsnvpos$count_div, allsnvpos$count_par+allsnvpos$count_div > 0)

# saveRDS(object = allsnvpos, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/allsnv_bialpos_wcounts_2021.RDS")
# saveRDS(object = allsnvsgr, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/allsnvsgr_2021.RDS")

tabdf <- data.frame(table(decode(allsnvpos$total), decode(allsnvpos$bial)), stringsAsFactors = F)
tabdf$Var2 <- as.logical(tabdf$Var2)
tabdf$nFreq <- tabdf$Freq / (tabdf[!tabdf$Var2, "Freq"] + tabdf[tabdf$Var2, "Freq"])

# View(as.data.frame(allsnvpos[which(allsnvpos$total == 4)], row.names = NULL))

p1 <- ggplot(data = tabdf, mapping = aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col(position = "dodge") + scale_y_log10() + theme_minimal()
p1

p1 <- ggplot(data = tabdf, mapping = aes(x = Var1, y = nFreq, fill = Var2)) + 
  geom_col() + theme_minimal()
p1



allsnvpos <- readRDS(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/allsnv_bialpos_wcounts_2021.RDS")
allsnvsgr <- readRDS(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/allsnvsgr_2021.RDS")

##### shrink objects to multiply hit loci only:
allsnvpos_multionly <- allsnvpos[which(allsnvpos$total > 1)]
allsnvsgr_multionly <- subsetByOverlaps(x = allsnvsgr, ranges = allsnvpos_multionly, type = "equal")
rm("allsnvpos", "allsnvsgr")

allsnvsgr_multionly$histology_abbreviation <- sumtab$histology_abbreviation[match(x = allsnvsgr_multionly$sampleid, table = sumtab$samplename)]
alldrivers <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
alldriversgr <- GRanges(seqnames = alldrivers$chr, ranges = IRanges(start = alldrivers$pos, width = nchar(alldrivers$ref)), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
alldriversgr$gene_name <- sapply(X = strsplit(alldrivers$gene_id, split = "::"), FUN = "[", 3)
alldriversgr$coding_noncoding <- alldrivers$coding_noncoding
alldriversgr <- unique(alldriversgr)

allsnvpos_multionly$gene_name <- ""
drivhits <- findOverlaps(query = allsnvpos_multionly, subject = alldriversgr, type = "equal")
mcols(allsnvpos_multionly)[queryHits(drivhits), "gene_name"] <- mcols(alldriversgr)[subjectHits(drivhits), "gene_name"]

gc35 <- import.gff3("/srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v35lift37.basic.annotation.gff3.gz", feature.type = "gene")
seqlevelsStyle(gc35) <- "ensembl"

allsnvpos_multionly$gene_name_gc35 <- ""
allsnvpos_multionly$gene_type <- ""
gc35hits <- findOverlaps(query = allsnvpos_multionly, subject = gc35)
mcols(allsnvpos_multionly)[queryHits(gc35hits), c("gene_name_gc35")] <- mcols(gc35)[subjectHits(gc35hits), c("gene_name")]

mcols(allsnvpos_multionly)[, c("gene_type", "gene_name_gc35_prox")] <- mcols(gc35)[nearest(x = allsnvpos_multionly, subject = gc35), c("gene_name", "gene_type")]

# tempdf <- as.data.frame(allsnvpos_multionly[head(order(allsnvpos_multionly$total, decreasing = T), n = 100)], row.names = NULL)
# tempdf <- as.data.frame(allsnvpos_multionly[allsnvpos_multionly$count_par > 0 | allsnvpos_multionly$count_div > 0 ], row.names = NULL)
# 
# allsnvsgr_multionly[which(start(allsnvsgr_multionly) == "87566156")]


# plotgr <- unique(c(allsnvpos_multionly[order(allsnvpos_multionly$count_par + allsnvpos_multionly$count_div, decreasing = T)][1:50],
#                    allsnvpos_multionly[which(allsnvpos_multionly$gene_name == ""), ][order(mcols(allsnvpos_multionly)[which(allsnvpos_multionly$gene_name == ""), "count"], decreasing = T)][1:50]))
# p1 <- ggplot(data = data.frame(count = allsnvpos_multionly$count, par = allsnvpos_multionly$count_par, div = allsnvpos_multionly$count_div))

### MOD to 52
# allsnvpos_multibial <- allsnvpos_multionly[allsnvpos_multionly$count_div + allsnvpos_multionly$count_par > 0] 
allsnvpos_multibial <- allsnvpos_multionly[order(allsnvpos_multionly$count_par + allsnvpos_multionly$count_div, allsnvpos_multionly$count, decreasing = T)][1:52]
allsnvsgr_multibial <- subsetByOverlaps(x = allsnvsgr_multionly, ranges = allsnvpos_multibial, type = "equal")
strand(allsnvsgr_multibial) <- ifelse(allsnvsgr_multibial$REF %in% c("C", "T"), "+", "-")

mcols(allsnvsgr_multibial)[, c("gene_type", "gene_name_gc35_prox")] <- mcols(gc35)[nearest(x = allsnvsgr_multibial, subject = gc35), c("gene_name", "gene_type")]
allsnvsgr_multibial$context <- getSeq(x = BSgenome.Hsapiens.1000genomes.hs37d5, names = resize(allsnvsgr_multibial, fix = "center", width = 15))

allsnvsgr_multibialdf <- as.data.frame(allsnvsgr_multibial, row.names = NULL)

tabout <- table(allsnvsgr_multibialdf$histology_abbreviation, paste0(allsnvsgr_multibialdf$seqnames, ":", allsnvsgr_multibialdf$start)) 
taboutdf <- as.data.frame(tabout[rowSums(tabout) > 1,])

#sort sites by type of mut (i.e. tumour type), then by freq of hits
# taboutdf$Var1 <- factor(taboutdf$Var1, levels = names(sort(rowSums(tabout), decreasing = T)))
taboutdf$Var1 <- factor(taboutdf$Var1, levels = c("Skin-Melanoma", "ColoRect-AdenoCA", "Uterus-AdenoCA", "Eso-AdenoCA", "Stomach-AdenoCA", "Liver-HCC", "Bladder-TCC", "Panc-AdenoCA", "Breast-AdenoCA", "Lung-SCC", "Thy-AdenoCA"))

orderv <- by(data = taboutdf, INDICES = taboutdf$Var1, function(x) {
  y <- x[order(x$Freq, decreasing = T),]
  y <- as.character(y[y$Freq > 0, "Var2"])
})
taboutdf$Var2 <- factor(taboutdf$Var2, levels = as.vector(unlist(orderv)[!duplicated(unlist(orderv))]))
# taboutdf$Var2 <- factor(taboutdf$Var2, levels = unique(taboutdf[order(taboutdf$Var1, taboutdf$Freq, decreasing = F), "Var2"]))




# 6:142706206 is sequence palindrome, seemingly hit by APOBEC
## additional plotting info, the numbers of par/div hits
parhits_tomerge <- subsetByOverlaps(x = parhits[, c("ref", "alt", "sampleid")], ranges = allsnvpos_multibial, type = "equal")
parhits_tomerge <- parhits_tomerge[rep(1:length(parhits_tomerge), rep(2, length(parhits_tomerge)))]
thirdhits_tomerge <- subsetByOverlaps(x = thirdhits[, c("REF", "alt1", "alt2", "sampleid")], ranges = allsnvpos_multibial, type = "equal")
colnames(mcols(thirdhits_tomerge)) <- c("ref", "alt", "alt2", "sampleid")
thirdhits_tomerge <- thirdhits_tomerge[rep(1:length(thirdhits_tomerge), rep(2, length(thirdhits_tomerge)))]
thirdhits_tomerge$alt <- ifelse(duplicated(paste0(seqnames(thirdhits_tomerge), ":", start(thirdhits_tomerge), "_", thirdhits_tomerge$sampleid)), thirdhits_tomerge$alt2, thirdhits_tomerge$alt)
thirdhits_tomerge <- thirdhits_tomerge[, c("ref", "alt", "sampleid")]
allbialhits_duplicated <- c(parhits_tomerge, thirdhits_tomerge)
allbialhits_duplicated$type <- factor(rep(c("par", "div"), c(length(parhits_tomerge), length(thirdhits_tomerge))))
allbialhits_duplicated$histology_abbreviation <- sumtab$histology_abbreviation[match(x = allbialhits_duplicated$sampleid, table = sumtab$samplename)]

tb1 <- by(data = allbialhits_duplicated$type, INDICES = paste0(allbialhits_duplicated$histology_abbreviation, "_", seqnames(allbialhits_duplicated), ":", start(allbialhits_duplicated)), FUN = function(x) c(table(x)))
tb1df <- data.frame(do.call(rbind, strsplit(names(tb1), split = "_")), do.call(rbind, tb1), row.names = NULL)
tb1df <- tb1df[tb1df$X2 %in% taboutdf$Var2, ] # this actually removes two false positives ... TODO: add further filtering by genotyping the normals at all thirdallele positions
# tb1df[!tb1df$X2 %in% taboutdf$Var2, ]

# taboutdf$Var2 <- factor(taboutdf$Var2, levels = c(levels(taboutdf$Var2), "1:17229152", "1:33475721"))


### add in mutational signatures
# for every sample, assign signature probabilities (i.e. each sig * exposures) and normalise to 1
monoatbial <- allsnvsgr_multibial[, c("REF", "ALT", "sampleid")]
colnames(mcols(monoatbial)) <- c("ref", "alt", "sampleid")
monoatbial$type <- "mono"
monoatbial$ref <- as.character(monoatbial$ref)
monoatbial$alt <- as.character(unlist(monoatbial$alt))
monoatbial$histology_abbreviation <- sumtab$histology_abbreviation[match(x = monoatbial$sampleid, table = sumtab$samplename)]
monoatbial <- monoatbial[!paste0(seqnames(monoatbial), ":", start(monoatbial), "_", monoatbial$histology_abbreviation) %in% 
                           paste0(seqnames(allbialhits_duplicated), ":", start(allbialhits_duplicated), "_", allbialhits_duplicated$histology_abbreviation), ]
monoatbial <- sort(c(monoatbial, allbialhits_duplicated))
strand(monoatbial) <- ifelse(monoatbial$ref %in% c("C", "T"), "+", "-")
monoatbial$context <- paste0(getSeq(x = BSgenome.Hsapiens.1000genomes.hs37d5, names = resize(monoatbial, fix = "center", width = 3)), ">", ifelse(strand(monoatbial) == "+", monoatbial$alt, reverseComplement(DNAStringSet(monoatbial$alt))))
monoatbial$largecontext <- getSeq(x = BSgenome.Hsapiens.1000genomes.hs37d5, names = resize(monoatbial, fix = "center", width = 15))

# prob for each sig, irrespective of exposures
sigs <- read.csv(file = "/srv/shared/vanloo/ICGC_signatures/20180322_release/sigProfiler_SBS_signatures.csv", as.is = T)
exposuresinsamples <- read.csv(file = "/srv/shared/vanloo/ICGC_signatures/20180322_release/PCAWG_sigProfiler_SBS_signatures_in_samples_waliqID.csv", as.is = T)
trinucfreq <- colSums(trinucleotideFrequency(x = getSeq(BSgenome.Hsapiens.1000genomes.hs37d5)[1:23]))
reftrinuc <- grep(pattern = "[ACGT][CT][ACGT]", x = names(trinucfreq), value = T)
trinucfreq <- trinucfreq[reftrinuc] + trinucfreq[as.character(reverseComplement(DNAStringSet(reftrinuc)))]

sigs$context <- sapply(X = strsplit(sigs$Mutation.Type, split = ">"), FUN = '[', 1)
sigs[, -c(1, 67)] <- sigs[, -c(1, 67)] * trinucfreq[match(x = sigs$context, table = names(trinucfreq))]
sigs[, -c(1, 67)] <- t(t(sigs[, -c(1, 67)]) / colSums(sigs[, -c(1, 67)]))

actweightedprobs <- sigs[match(x = monoatbial$context, table = sigs$Mutation.Type), -c(1, 67)] * exposuresinsamples[match(x = monoatbial$sampleid, table = exposuresinsamples$Sample.Name), -c(1:3)]
actweightedprobs <- actweightedprobs / rowSums(actweightedprobs)
# rowSums(actweightedprobs)
# (sigspercontext * exposurespersample)[1:10, 1:10]

actweightedprobsperpos <- data.frame(do.call(rbind, by(data = actweightedprobs, INDICES = paste0(seqnames(monoatbial), ":", start(monoatbial)), FUN = function(x) colSums(x)/sum(x))))
actweightedprobsperpos["15:74262402",]
actweightedprobsperpos$posid <- rownames(actweightedprobsperpos)
actweightedprobsperpos$largecontext <- as.character(monoatbial$largecontext[match(x = actweightedprobsperpos$posid, table = paste0(seqnames(monoatbial), ":", start(monoatbial)))])
sort(colSums(actweightedprobsperpos[,-c(66,67)]))

p1 <- ggplot(data = taboutdf, mapping = aes(x = Var1, y = Var2)) + geom_raster(mapping = aes(fill = Freq)) + scale_fill_gradient(low = "#ffffff", high = "#d73027")
p1 <- p1 + geom_text(data = tb1df, mapping = aes(x = X1, y = X2, label = paste0(par/2, "/", div/2)))
p1 <- p1 + theme_minimal() + scale_y_discrete(limits = rev(levels(taboutdf$Var2))) + theme(axis.text.x = element_text(angle = 90))
p1
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Rheinbay_type_recurrently_biallelic_sites_heatmap_50_2021.pdf"), plot = p1, width = 5, height = 10, useDingbats=FALSE)

p1 <- ggplot(data = taboutdf[taboutdf$Var1 == "Skin-Melanoma", ], mapping = aes(x = Var1, y = Var2, label = Freq)) + geom_ba() + geom_raster(mapping = aes(fill = Freq)) + scale_fill_gradient(low = "#f7f7f7", high = "#67001f") + theme_minimal()
p1

library(reshape2)
actweightedprobsperposlong <- melt(data = actweightedprobsperpos, id.vars = c("posid", "largecontext"))
actweightedprobsperposlong$variable <- as.character(actweightedprobsperposlong$variable)
actweightedprobsperposlong$variable[actweightedprobsperposlong$variable %in% paste0("SBS", c(1,5,40))] <-  "Age"
actweightedprobsperposlong$variable[actweightedprobsperposlong$variable %in% paste0("SBS", c(2,13))] <-  "APOBEC"
actweightedprobsperposlong$variable[actweightedprobsperposlong$variable %in% paste0("SBS", c(10), c("a","b"))] <-  "POLE"
actweightedprobsperposlong$variable[actweightedprobsperposlong$variable %in% paste0("SBS", 17, c("a","b"))] <-  "Sig17"
actweightedprobsperposlong$variable[actweightedprobsperposlong$variable %in% paste0("SBS", 7, c("a", "b", "c", "d"))] <-  "UV"
actweightedprobsperposlong$variable[actweightedprobsperposlong$variable %in% paste0("SBS", c(6,14,15,20,21,26,44))] <-  "MSI"
actweightedprobsperposlong$variable[grep(pattern = "SBS", x = actweightedprobsperposlong$variable)] <-  "Other"
actweightedprobsperposlong$variable <- factor(x = actweightedprobsperposlong$variable, levels = c("Age", "APOBEC", "MSI", "POLE", "Sig17", "UV", "Other"))
p2 <- ggplot(data = actweightedprobsperposlong, mapping = aes(x = posid)) + geom_col(mapping = aes(y = value, fill = variable)) + scale_x_discrete(limits = rev(levels(taboutdf$Var2))) + coord_flip(clip = "off", ylim = c(0, 3)) + theme_minimal() + theme(panel.grid = element_blank(), axis.text.y = element_blank())
p2 <- p2 + geom_text(mapping = aes(label = largecontext, y = 2))
p2 <- p2 + scale_fill_manual(values = c(Age = "#8dd3c7", MSI= "#ffffb3", POLE = "#bebada", APOBEC = "#fb8072", UV = "#80b1d3", Sig17 = "#fdb462", Other = "#D3D3D3"))
p2
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/Rheinbay_type_recurrently_biallelic_sites_marginalbars_50_2021.pdf"), plot = p2, width = 4.5, height = 10, useDingbats=FALSE)

allsnvpos_multionly$gene_name_gc35
# subset all snvs to positions with a few hits
# compute tables like sup fig for this.



allsnvpos[allsnvpos$total %in% c(20, 29)]
allsnvpos[which(allsnvpos$total == 1 & allsnvpos$bial)]

allsnvposdf <- data.frame(mcols(allsnvpos)[, c("Variant_Classification", "count")])

allpos <- c(granges(allsnvsgr), granges(parhits), granges(thirdhits))
multihitpos <- allpos[overlapsAny(query = allpos, type = "equal", drop.self = T, drop.redundant = T)]

multihitpos <- allsnvpos[allsnvpos$count > 1]
multihitpossnvs <- subsetByOverlaps(x = allsnvsgr, ranges = multihitpos, type = "equal")

tabledf <- table(allsnvpos$count)

p1 <- ggplot(data = allsnvposdf, mapping = aes(x = count)) + geom_bar() + scale_y_log10() + labs(x = )
p1

thirdhitsforaddition <- granges(thirdhits)
mcols(thirdhitsforaddition) <- DataFrame(alt = thirdhits$alt1,
                                         trinuc = getSeq(resize(thirdhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5),
                                         ref = thirdhits$REF)
# thirdhitsforaddition <- c(granges(thirdhits),granges(thirdhits))
# mcols(thirdhitsforaddition) <- DataFrame(alt = c(thirdhits$alt1,thirdhits$alt2),
#                                          trinuc = rep(getSeq(resize(thirdhits, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5), 2),
#                                          ref = rep(thirdhits$REF, 2))
flipthirdidx <- which(thirdhitsforaddition$ref %in% c("A", "G"))
thirdhitsforaddition$trinuc[flipthirdidx] <- reverseComplement(thirdhitsforaddition$trinuc[flipthirdidx])
thirdhitsforaddition$alt <- DNAStringSet(thirdhitsforaddition$alt)
thirdhitsforaddition$alt[flipthirdidx] <- reverseComplement(thirdhitsforaddition$alt[flipthirdidx])


driverhits <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
driverhits <- driverhits[driverhits$ref %in% c("A", "C", "G", "T") & driverhits$alt %in% c("A", "C", "G", "T"), c("chr", "pos")]
driverhits <- sort(unique(GRanges(seqnames = driverhits$chr, ranges = IRanges(start = driverhits$pos, width = 1))))

# sum(sapply(sampleids, FUN = function(x, snvsmnvdir) file.exists(file.path(snvsmnvdir, pattern = paste0(x, ".consensus.20160830.somatic.snv_mnv.vcf.gz"))), snvsmnvdir = SNVMNVINDELDIR))
# 2583/2658 samples

get_snvs <- function(sampleid, snvsmnvdir) {
  print(paste0("Loading snvs for sample ", sampleid))
  snv_mnvfile <- file.path(snvsmnvdir, pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz"))
  
  if (!file.exists(snv_mnvfile)) {
    snvs <- GRanges()
    snvs$alt <- character()
  } else {
    snvs <- rowRanges(readVcf(file = snv_mnvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
    snvs <- snvs[lengths(mcols(snvs)$ALT) == 1 & lengths(mcols(snvs)$REF) == 1 & seqnames(snvs) != "Y"]
    mcols(snvs) <- DataFrame(alt = as.character(unlist(snvs$ALT)))
  }
  
  snvs$trinuc <- getSeq(resize(snvs, fix = "center", width = 3), x = BSgenome.Hsapiens.1000genomes.hs37d5)
  snvs$ref <- subseq(snvs$trinuc, 2,2)
  flipidx <- which(snvs$ref %in% c("A", "G"))
  snvs$trinuc[flipidx] <- reverseComplement(snvs$trinuc[flipidx])
  snvs$alt <- DNAStringSet(snvs$alt)
  snvs$alt[flipidx] <- reverseComplement(snvs$alt[flipidx])
  snvs$sampleid <- rep(sampleid, length(snvs))
  # snvs <- snvs[which(snvs$trinuc == "TCT")]
  
  return(snvs)
}

### remove all driver muts to reduce influence of selection
# debug(get_snvs)
allsnvs <- unlist(GRangesList(mclapply(X = sampleids, FUN = get_snvs, snvsmnvdir = SNVMNVINDELDIR, mc.preschedule = T, mc.cores = 16)))
# allsnvs <- c(allsnvs, thirdhitsforaddition)#[which(thirdhitsforaddition$trinuc == "TCT")])
# allsnvs$trinucmut <- factor(paste0(allsnvs$trinuc, ">", allsnvs$alt))
# allsnvs <- subsetByOverlaps(x = allsnvs, ranges = driverhits, invert = T)
# 