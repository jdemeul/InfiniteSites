## check for each signature probs of forward, back (same copy) and parallel, third allele (other copy)

library(ggplot2)
library(Biostrings)

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")

bases <- generate_bases_types_trinuc()
bases$trinucleotides_mutations_gt2 <- paste0(substr(bases$trinucleotides_mutations_gt, 5,7),substr(bases$trinucleotides_mutations_gt, 2,3))
bases$trinucleotides_mutations_gt2 <- factor(bases$trinucleotides_mutations_gt2, levels = bases$trinucleotides_mutations_gt2)

PCAWGSIGFILE <- "/srv/shared/vanloo/ICGC_signatures/20180322_release/sigProfiler_SBS_signatures.csv"

sigs <- read.csv(file = PCAWGSIGFILE, as.is = T)

rownames(sigs) <- sigs$Mutation.Type

longdf <- expand.grid(sigs$Mutation.Type, sigs$Mutation.Type, stringsAsFactors = F)
colnames(longdf) <- c("mut1", "mut2")


longdf[, c("reftri1", "alt1")] <- do.call(rbind, strsplit(longdf$mut1, split = ">"))
longdf$ref1 <- substr(x = longdf$mut1, 2, 2)
longdf$alttri1 <- paste0(substr(x = longdf$mut1, 1, 1), longdf$alt1, substr(x = longdf$mut1, 3, 3))
longdf$alttri1 <- ifelse(longdf$alt1 %in% c("A", "G"), as.character(reverseComplement(DNAStringSet(longdf$alttri1))), longdf$alttri1)


longdf[, c("reftri2", "alt2")] <- do.call(rbind, strsplit(longdf$mut2, split = ">"))
longdf$ref2 <- substr(x = longdf$mut2, 2, 2)
longdf$alttri2 <- paste0(substr(x = longdf$mut2, 1, 1), longdf$alt2, substr(x = longdf$mut2, 3, 3))
longdf$alttri2 <- ifelse(longdf$alt1 %in% c("A", "G"), as.character(reverseComplement(DNAStringSet(longdf$alttri2))), longdf$alttri2)


longdf$mut1 <- factor(longdf$mut1, levels = levels(bases$trinucleotides_mutations_gt2))
longdf$mut2 <- factor(longdf$mut2, levels = levels(bases$trinucleotides_mutations_gt2))

#empty/impossible unless stated otherwise
longdf$type <- "empty"

# parallel defined as: 2x same reftri + alt
# third allele: 2x same reftri + different alt

# forward: alttri = reftri & alt2 != ref1
# back: alttri = reftri & alt2 = ref1

longdf$type <- ifelse(longdf$reftri1 == longdf$reftri2,
                      ifelse(longdf$alt1 == longdf$alt2, "parallel", "third_allele"),
                      ifelse(longdf$alttri1 == longdf$reftri2, ifelse(longdf$alttri2 == longdf$reftri1, "back", "forward"), longdf$type))
longdf$type <- factor(longdf$type, levels = c("back", "forward", "parallel", "third_allele", "empty"))
longdf <- longdf[!(longdf$type == "third_allele" & as.integer(longdf$mut1) > as.integer(longdf$mut2)), ]



OUTDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isaviolationtypes_persig/"

signats <- grep(pattern = "^SBS", x = colnames(sigs), value = T)

for (signat in signats) {
  longdf$prob <- exp(log(sigs[as.character(longdf$mut1), signat])+log(sigs[as.character(longdf$mut2), signat]))
  longdf[longdf$type %in% c("parallel", "third_allele"), "prob"] <- 2*longdf[longdf$type %in% c("parallel", "third_allele"), "prob"]

  longdf[longdf$type == "empty", "prob"] <- 0
  longdf$prob <- longdf$prob/sum(longdf$prob)

  
  by(data = longdf$prob, INDICES = longdf$type, FUN = sum)
  
  p1 <- ggplot(data = longdf[longdf$type != "empty", ], mapping = aes(x = mut1, y = mut2)) + geom_tile(mapping = aes(fill = type, alpha = ifelse(prob >= .001, 1, .1)), show.legend = F)
  p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + scale_fill_brewer(type = "qual", palette = "Set1") #+ scale_alpha_discrete(values = c('TRUE' = 1, 'FALSE' = .1))
  # p1
  
  p2 <- ggplot(data = longdf[longdf$type != "empty", ], mapping = aes(x = type, fill = type)) + geom_bar(mapping = aes(weight = prob)) 
  p2 <- p2 + theme_minimal() + scale_fill_brewer(type = "qual", palette = "Set1")
  # p2
  
  ggsave(filename = file.path(OUTDIR, paste0(signat, "_heatmap_p01.pdf")), plot = p1, width = 12, height = 12)
  ggsave(filename = file.path(OUTDIR, paste0(signat, "_typeprobs.pdf")), plot = p2, width = 4, height = 4)
}
