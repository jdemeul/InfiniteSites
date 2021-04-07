### check expression in MELA samples with and without RPL18A biallelic hits etc

expr <- read.delim(file = "/srv/shared/vanloo/ICGC-gene_expression/joint_fpkm_uq.tsv.gz", nrows = 10)
sumtab <- read.delim(file = "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv", as.is = T)
sumtab <- sumtab[sumtab$dcc_project_code == "SKCM-US" & sumtab$wgs_exclusion_white_gray == "Whitelist" & sumtab$tumor_rna_seq_aliquot_id != "", ]
infsitesmut <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/TopBialHitSite.txt", as.is = T)

colv <- rep("NULL", length(colnames(expr)))
colv[1] <- "character"
colv[which(colnames(expr) %in% paste0("X", gsub(pattern = "-", replacement = ".", x = sumtab$tumor_rna_seq_aliquot_id)))] <- "numeric"

expr <- read.delim(file = "/srv/shared/vanloo/ICGC-gene_expression/joint_fpkm_uq.tsv.gz", colClasses = colv)

colnames(expr) <- c("feature", sumtab$tumor_wgs_aliquot_id[match(x = colnames(expr)[-1], table = paste0("X", gsub(pattern = "-", replacement = ".", x = sumtab$tumor_rna_seq_aliquot_id)))])

expr$feature <- sapply(X = strsplit(x = expr$feature, split = "::"), FUN = "[", 2)

sort(unlist(expr[grep(pattern = "RPL18A::", x = expr$feature), -1]))

library(reshape2)

exprmelt <- melt(data = expr, id.vars = "feature", variable.name = "sampleid", value.name = "fpkm")
head(exprmelt)

exprmelt$ishit <- infsitesmut$type[match(x = exprmelt$sampleid, table = infsitesmut$sampleid)]
exprmelt$ishit[is.na(exprmelt$ishit)] <- "none"

library(ggplot2)

p1 <- ggplot(data = exprmelt[grep(pattern = "RPL18", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1

p1 <- ggplot(data = exprmelt[grep(pattern = "KDM5A", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1

p1 <- ggplot(data = exprmelt[grep(pattern = "NFAT|ET[SV]|ELK", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1

p1 <- ggplot(data = exprmelt[grep(pattern = "RPS20", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1

p1 <- ggplot(data = exprmelt[grep(pattern = "KDM5A", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1

p1 <- ggplot(data = exprmelt[grep(pattern = "RP5-1125", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1

p1 <- ggplot(data = exprmelt[grep(pattern = "IQGAP1", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1

p1 <- ggplot(data = exprmelt[grep(pattern = "GSPT1", x = exprmelt$feature), ]) + geom_jitter(width = 0.2, height = 0, mapping = aes(x = feature, y = log2(fpkm), colour = ishit))
p1
