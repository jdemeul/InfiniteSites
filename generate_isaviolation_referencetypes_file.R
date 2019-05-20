### generate template file for ISA violations

library(Biostrings)

generate_bases_types_trinuc <- function() {
  # generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  types <- c("A", "G", "T", "A", "C", "G")
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(trinucleotides, ">", rep(types, rep(16,6))))
  return(list(bases = bases, types = types, trinucleotides = trinucleotides,
              trinucleotides_mutations = trinucleotides_mutations))
}


bases <- generate_bases_types_trinuc()

longdf <- expand.grid(bases$trinucleotides_mutations, bases$trinucleotides_mutations, stringsAsFactors = F, KEEP.OUT.ATTRS = F)
# all(levels(longdf$Var1) == levels(longdf$Var2))
colnames(longdf) <- c("mut1", "mut2")


longdf[, c("reftri1", "alt1")] <- do.call(rbind, strsplit(longdf$mut1, split = ">"))
longdf$ref1 <- substr(x = longdf$mut1, 2, 2)
longdf$alttri1 <- paste0(substr(x = longdf$mut1, 1, 1), longdf$alt1, substr(x = longdf$mut1, 3, 3))
longdf$alttri1 <- ifelse(longdf$alt1 %in% c("A", "G"), as.character(reverseComplement(DNAStringSet(longdf$alttri1))), longdf$alttri1)


longdf[, c("reftri2", "alt2")] <- do.call(rbind, strsplit(longdf$mut2, split = ">"))
longdf$ref2 <- substr(x = longdf$mut2, 2, 2)
longdf$alttri2 <- paste0(substr(x = longdf$mut2, 1, 1), longdf$alt2, substr(x = longdf$mut2, 3, 3))
longdf$alttri2 <- ifelse(longdf$alt1 %in% c("A", "G"), as.character(reverseComplement(DNAStringSet(longdf$alttri2))), longdf$alttri2)

longdf$mut1 <- factor(longdf$mut1, levels = bases$trinucleotides_mutations)
longdf$mut2 <- factor(longdf$mut2, levels = bases$trinucleotides_mutations)

#empty/impossible unless stated otherwise
longdf$idx <- 1:nrow(longdf)
longdf$type <- NA

# parallel defined as: 2x same reftri + alt
# third allele: 2x same reftri + different alt

# forward: alttri = reftri & alt2 != ref1
# back: alttri = reftri & alt2 = ref1

longdf$type <- ifelse(longdf$reftri1 == longdf$reftri2,
                      ifelse(longdf$alt1 == longdf$alt2, "parallel", "third_allele"),
                      ifelse(longdf$alttri1 == longdf$reftri2, ifelse(longdf$alttri2 == longdf$reftri1, "back", "forward"), longdf$type))
longdf <- longdf[!((longdf$type == "third_allele" & as.integer(longdf$mut1) > as.integer(longdf$mut2)) | is.na(longdf$type)), ]
longdf$type <- factor(longdf$type, levels = c("back", "forward", "parallel", "third_allele"))

write.table(x = longdf[, c("idx", "mut1", "mut2", "type")], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isavioltypes_reffile.txt", quote = F, col.names = T, sep = "\t", row.names = F)


