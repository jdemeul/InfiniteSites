

snv_phasing = function(germ_loci_file = NA, som_vcf_file = NA, outfile, bam_file, bai_file, max_distance, minMapQ = 20, minBaseQ = 20, max_mm_read = 2, max_mm_pair = 3, chrominfo = genome_seqinfo, ncores = 12, preschedule = F, chr_subset = NULL) {
  
  count_data <- data.frame(chr = character(), pos1 = integer(), ref1 = character(), alt1 = character(),
                           pos2 = integer(), ref2 = character(), alt2 = character(),
                           Num_ref_ref = integer(), Num_alt_alt = integer(), Num_alt_ref = integer(), Num_ref_alt = integer(),
                           phasing = character(), stringsAsFactors = F)
  snvs <- data.frame(chr = character(), pos = integer(), ref = character(), alt = character(), type = character(), stringsAsFactors = F)
  
  if (!is.na(germ_loci_file) && file.exists(germ_loci_file)) {
    # gl_snvs <- read.delim(germ_loci_file, header = F, row.names = NULL, stringsAsFactors = F, col.names = c("chr", "pos", "ref", "alt"))
    gl_snvs <- read_tsv(file = germ_loci_file, col_names = T, col_types = "cicc")
    gl_snvs$type <- "SNP"
    snvs <- rbind(snvs, gl_snvs)
    rm(gl_snvs)
  }
  if (!is.na(som_vcf_file) && file.exists(som_vcf_file)) {
    # som_snvs <- read.delim(som_loci_file, header = F, row.names = NULL, stringsAsFactors = F, col.names = c("chr", "pos", "ref", "alt"))
    # som_snvs <- read_tsv(file = som_loci_file, col_names = c("chr", "pos", "ref", "alt"), col_types = "cicc")
    som_vcf <- SummarizedExperiment::rowRanges(VariantAnnotation::readVcf(file = som_vcf_file), genome = chrominfo)
    som_snvs <- data.frame(chr = as.character(seqnames(som_vcf)), pos = start(som_vcf),
                         ref = as.character(mcols(som_vcf)$REF), alt = as.character(unlist(mcols(som_vcf)$ALT)), type = "SNV")
    snvs <- rbind(snvs, som_snvs)
    rm(som_snvs, som_vcf)
  }
  
  if (nrow(snvs) > 0) {
    # fai <- parseFai(fai_file)
    # ign <- parseIgnore(ign_file)
    # fai <- fai[!(fai$chromosome %in% ign$chromosome), ]
    
    # chrominfo <- Seqinfo(seqnames = fai$chromosome, seqlengths = fai$length, isCircular = rep(F, nrow(fai)), genome = "GRCh37")
    
    snvs_gr <- GRanges(seqnames = snvs$chr, ranges = IRanges(start = snvs$pos, end = snvs$pos), mcols = snvs[, -c(1,2)], seqinfo = chrominfo)
    # if (!is.null(keep_chroms)) {
    #   snvs_gr <- keepSeqlevels(x = snvs_gr, value = keep_chroms, pruning.mode = "coarse")
    # }
    snvs_gr <- sort(snvs_gr)
    # for testing purposes
    if (!is.null(chr_subset)) {
        snvs_gr <- keepSeqlevels(x = snvs_gr, value = chr_subset, pruning.mode = "coarse")
    }
    
    ## start parallelising here, do the following per chromosome
    # if (parallel) {
      count_data <- mclapply(X = split(snvs_gr, f = seqnames(snvs_gr)), FUN = snv_phasing_chr, minBaseQ = minBaseQ, minMapQ = minMapQ, max_mm_read = max_mm_read, max_mm_pair = max_mm_pair, bam_file = bam_file, bai_file = bai_file, mc.preschedule = preschedule, mc.cores = ncores)
    # } else {
    #   count_data <- lapply(X = split(snvs_gr, f = seqnames(snvs_gr)), FUN = snv_phasing_chr, minBaseQ = minBaseQ, minMapQ = minMapQ, max_mm_read = max_mm_read, max_mm_pair = max_mm_pair, bam_file = bam_file, bai_file = bai_file)
    # }
    count_data <- do.call(rbind, count_data)
    
  }
  write.table(x = count_data, file = outfile, sep = "\t", quote = F, row.names = F)
  return(NULL)
}



snv_phasing_chr <- function(snvs_gr, minBaseQ, minMapQ, max_mm_read, max_mm_pair, bam_file, bai_file) {
  count_data <- data.frame(chr = character(), pos1 = integer(), ref1 = character(), alt1 = character(), type1 = character(),
                           pos2 = integer(), ref2 = character(), alt2 = character(), type2 = character(),
                           Num_ref_ref = integer(), Num_alt_alt = integer(), Num_alt_ref = integer(), Num_ref_alt = integer(),
                           phasing = character(), stringsAsFactors = F)
  
  if(length(snvs_gr) > 0) {
    proximityHits <- findOverlaps(query = snvs_gr, maxgap = max_distance, drop.self = T, drop.redundant = T)
    snvs_gr1 <- snvs_gr[queryHits(proximityHits)]
    snvs_gr2 <- snvs_gr[subjectHits(proximityHits)]
    
    if (length(snvs_gr1) > 0) {
      # ref_ref = paste(mcols(snvs_gr1)$mcols.ref, mcols(snvs_gr2)$mcols.ref, sep = "_")
      # alt_ref = paste(mcols(snvs_gr1)$mcols.alt, mcols(snvs_gr2)$mcols.ref, sep = "_")
      # ref_alt = paste(mcols(snvs_gr1)$mcols.ref, mcols(snvs_gr2)$mcols.alt, sep = "_")
      # alt_alt = paste(mcols(snvs_gr1)$mcols.alt, mcols(snvs_gr2)$mcols.alt, sep = "_")
      # 
      # avoid scanning the same region twice
      rep_idxs <- match(start(snvs_gr1), unique(start(snvs_gr1)))
      
      # flag = scanBamFlag(isPaired = T, hasUnmappedMate = F, isDuplicate = F, isUnmappedQuery = F) #, isProperPair=T
      flag = scanBamFlag(isNotPassingQualityControls = F, isDuplicate = F, isProperPair = T)
      param <- ScanBamParam(which = snvs_gr1[!duplicated(start(snvs_gr1))],
                            what = c("groupid", "pos", "mpos", "qwidth", "seq", "qual"),
                            flag = flag, mapqFilter = minMapQ, tag = c("MD"), simpleCigar = T) #, mapqFilter=minMapQ # simpleCIGAR: allow only reads with M in CIGAR string, no clipping, indels or whatever
      bamfile <- BamFile(file = bam_file, index = bai_file, asMates = T)
      bam <- scanBam(bamfile, param = param)
      count_data <- t(mapply(FUN = get_allele_combination_counts, bam = bam[rep_idxs], snv1pos = start(snvs_gr1), snv2pos = start(snvs_gr2), 
                        snv1type = mcols(snvs_gr1)$mcols.type, snv2type = mcols(snvs_gr2)$mcols.type, snv1ref = mcols(snvs_gr1)$mcols.ref,
                        snv1alt = mcols(snvs_gr1)$mcols.alt, snv2ref = mcols(snvs_gr2)$mcols.ref, snv2alt = mcols(snvs_gr2)$mcols.alt,
                        MoreArgs = list(minBaseQ = minBaseQ, max_mm_read = max_mm_read, max_mm_pair = max_mm_pair), SIMPLIFY = T))
      # alleles <- mapply(FUN = get_allele_combination_counts_gr, bam = bam, snv1pos = start(snvs_gr1), snv2pos = start(snvs_gr2), MoreArgs = list(minBaseQ = minBaseQ), SIMPLIFY = F)
      # alleles = get_allele_combination_counts_(bam, dat_pair)
      
      # count_data <- t(mapply(FUN = function(allele_combs, ref_ref, alt_ref, ref_alt, alt_alt) table(factor(paste(allele_combs$snv1, allele_combs$snv2, sep = "_"),
      #                                                                                                      levels = c(ref_ref, alt_alt, alt_ref, ref_alt),
      #                                                                                                      labels = c("Num_ref_ref", "Num_alt_alt", "Num_alt_ref", "Num_ref_alt"))),
      #                        allele_combs = alleles, ref_ref = ref_ref, alt_ref = alt_ref, ref_alt = ref_alt, alt_alt = alt_alt))
      
      count_data <- data.frame(chr = seqnames(snvs_gr1), pos1 = start(snvs_gr1), ref1 = mcols(snvs_gr1)$mcols.ref, alt1 = mcols(snvs_gr1)$mcols.alt, type1 = mcols(snvs_gr1)$mcols.type,
                               pos2 = start(snvs_gr2), ref2 = mcols(snvs_gr2)$mcols.ref, alt2 = mcols(snvs_gr2)$mcols.alt, type2 = mcols(snvs_gr2)$mcols.type,
                               count_data, stringsAsFactors = F, row.names = NULL)
      
      count_data <- count_data[rowSums(count_data[, c("Num_ref_ref", "Num_alt_alt", "Num_alt_ref", "Num_ref_alt")]) > 0, ]
      
      count_data$phasing <- annotate_phasing(count_data)
    }
  }
  return(count_data)
}



annotate_phasing <- function(countdata) {
  
  phasing <- ifelse(countdata$type1 == "SNP",
                    ifelse(countdata$type2 == "SNP",
                           ifelse(countdata$Num_alt_alt + countdata$Num_ref_ref > 0, # SNP-SNP AA+RR>0
                                  ifelse(countdata$Num_ref_alt + countdata$Num_alt_ref == 0, "phased", "converted"),
                                  ifelse(countdata$Num_ref_alt + countdata$Num_alt_ref > 0, "antiphased", "uninformative")),
                           ifelse(countdata$Num_alt_alt > 0,  # SNP-SNV AA>0
                                  ifelse(countdata$Num_ref_alt == 0, "phased", "converted/infsites"),
                                  ifelse(countdata$Num_ref_alt > 0, "antiphased", "uninformative"))),
                    ifelse(countdata$type2 == "SNV",
                           ifelse(countdata$Num_alt_alt > 0,
                                  ifelse(countdata$Num_ref_alt > 0,
                                         ifelse(countdata$Num_alt_ref == 0, "phased_subclone-clone", "converted/infsites"), # SNV-SNV AA>0, RA>0
                                         ifelse(countdata$Num_alt_ref > 0, "phased_clone-subclone", "phased_na")), #SNV-SNV AA>0, RA=0
                                  ifelse(countdata$Num_alt_ref + countdata$Num_ref_alt > 0, "antiphased_copies/subclones", "uninformative")), # SNV-SNV AA=0, RA+AR>0, informative if in CN = 1 region
                           ifelse(countdata$Num_alt_alt > 0, # SNV-SNP AA>0
                                  ifelse(countdata$Num_alt_ref == 0, "phased", "converted/infsites"),
                                  ifelse(countdata$Num_alt_ref > 0, "antiphased", "uninformative"))))
  
  return(phasing)
}



get_allele_combination_counts = function(bam, snv1pos, snv2pos, snv1type, snv2type, snv1ref, snv1alt, snv2ref, snv2alt, minBaseQ, max_mm_read, max_mm_pair) {
  # No two reads in pair when groupid is not dupicated - possibly removed due to mapping quality constraints
  # Likewise, sometimes a single read is present twice (3x same groupid)
  groupidcounts <- table(bam$groupid)
  proper_groupids <- as.integer(names(groupidcounts)[groupidcounts == 2])
  # exclude reads with indels or hard (e.g. SV)/soft clipped bases (now taken into account by simpleCIGAR = T in ScanBamParam)
  # proper_groupids <- setdiff(proper_groupids, bam$groupid[grep(pattern = "[IDSH]", x = bam$cigar)])
  
  if (length(proper_groupids) == 0) {
    return(c(Num_ref_ref = 0, Num_alt_alt = 0, Num_alt_ref = 0, Num_ref_alt = 0))
  }
  
  alleles <- as.data.frame(matrix(data = NA, nrow = length(proper_groupids), ncol = 14), stringsAsFactors = F)
  colnames(alleles) <- c("groupids", "j", "rel_pos_snv1", "read_snv1", "rel_pos_snv2", "read_snv2",
                         "read_1", "read_2", "qual_base1", "qual_base2", "base_snv1", "base_snv2", "mm_read1", "mm_read2")
  # alleles <- data.frame(groupids = proper_groupids, stringsAsFactors = F)
  alleles$groupids <- proper_groupids
  alleles$j <- match(alleles$groupids, bam$groupid)
  # 
  # alleles[ , c("rel_pos_snv1", "read_snv1", "rel_pos_snv2", "read_snv2", "read_1", "read_2", "qual_base1", "qual_base2", "base_snv1", "base_snv2", "mm_read1", "mm_read2")] <- NA
  
  # look at SNV 1
  rel_pos_read1_snv1 = (snv1pos-bam$pos[alleles$j]) + 1
  rel_pos_read2_snv1 = (snv1pos-bam$mpos[alleles$j]) + 1
  
  alleles$read_1 <- ifelse(rel_pos_read1_snv1 < bam$qwidth[alleles$j] & rel_pos_read1_snv1 > 0, T,
                           ifelse(rel_pos_read2_snv1 < bam$qwidth[alleles$j] & rel_pos_read2_snv1 > 0, F, NA))
  alleles$rel_pos_snv1 <- ifelse(alleles$read_1, rel_pos_read1_snv1, rel_pos_read2_snv1)
  alleles$read_snv1 <- ifelse(alleles$read_1, alleles$j, alleles$j + 1)
  
  # look at SNV 2  
  rel_pos_read1_snv2 = (snv2pos-bam$pos[alleles$j]) + 1
  rel_pos_read2_snv2 = (snv2pos-bam$mpos[alleles$j]) + 1
  
  alleles$read_2 <- ifelse(rel_pos_read1_snv2 < bam$qwidth[alleles$j] & rel_pos_read1_snv2 > 0, T,
                           ifelse(rel_pos_read2_snv2 < bam$qwidth[alleles$j] & rel_pos_read2_snv2 > 0, F, NA))
  alleles$rel_pos_snv2 <- ifelse(alleles$read_2, rel_pos_read1_snv2, rel_pos_read2_snv2)
  alleles$read_snv2 <- ifelse(alleles$read_2, alleles$j, alleles$j + 1)
  
  # get number of mismatches in reads
  alleles[, c("mm_read1", "mm_read2")] <- matrix(data = lengths(regmatches(x = bam$tag$MD, m = gregexpr("[ACGT]+", bam$tag$MD)))[bam$groupid %in% proper_groupids], ncol = 2, byrow = T)
  
  alleles <- alleles[ !(is.na(alleles$rel_pos_snv1) | is.na(alleles$rel_pos_snv2)), ]
  
  if (nrow(alleles) > 0) {
    alleles$qual_base1 <- mapply(FUN = function(qs, relpos) qs[relpos], qs = as(bam$qual[alleles$read_snv1], "IntegerList"), relpos = alleles$rel_pos_snv1)
    alleles$base_snv1 <- mapply(FUN = function(dnas, relpos) substr(as.character(dnas), relpos, relpos), dnas = bam$seq[alleles$read_snv1], relpos = alleles$rel_pos_snv1)
    
    alleles$qual_base2 <- mapply(FUN = function(qs, relpos) as.integer(qs)[relpos], qs = as(bam$qual[alleles$read_snv2], "IntegerList"), relpos = alleles$rel_pos_snv2)
    alleles$base_snv2 <- mapply(FUN = function(dnas, relpos) substr(as.character(dnas), relpos, relpos), dnas = bam$seq[alleles$read_snv2], relpos =alleles$ rel_pos_snv2)
    
    if (snv1type == "SNP") {
      alleles$mm_read1 <- ifelse(alleles$base_snv1 == snv1alt,
                                 ifelse(alleles$j == alleles$read_snv1, alleles$mm_read1 - 1, alleles$mm_read1),
                                 alleles$mm_read1)
      alleles$mm_read2 <- ifelse(alleles$base_snv1 == snv1alt,
                                 ifelse(alleles$j != alleles$read_snv1, alleles$mm_read2 - 1, alleles$mm_read2),
                                 alleles$mm_read2)
    }
    
    if (snv2type == "SNP") {
      alleles$mm_read1 <- ifelse(alleles$base_snv2 == snv2alt,
                                 ifelse(alleles$j == alleles$read_snv2, alleles$mm_read1 - 1, alleles$mm_read1),
                                 alleles$mm_read1)
      alleles$mm_read2 <- ifelse(alleles$base_snv2 == snv2alt,
                                 ifelse(alleles$j != alleles$read_snv2, alleles$mm_read2 - 1, alleles$mm_read2),
                                 alleles$mm_read2)
    }
    
    alleles <- alleles[alleles$qual_base1 > minBaseQ & alleles$qual_base2 > minBaseQ &
                         alleles$mm_read1 + alleles$mm_read2 <= max_mm_pair &
                         alleles$mm_read1 <= max_mm_read & alleles$mm_read2 <= max_mm_read, ]
  } else {
    return(c(Num_ref_ref = 0, Num_alt_alt = 0, Num_alt_ref = 0, Num_ref_alt = 0))
  }
  
  if (nrow(alleles) == 0) {
    return(c(Num_ref_ref = 0, Num_alt_alt = 0, Num_alt_ref = 0, Num_ref_alt = 0))
  } else {
    count_data <- c(table(factor(paste0(alleles$base_snv1, alleles$base_snv2),
                                 levels = c(paste0(snv1ref, snv2ref), paste0(snv1alt, snv2alt), paste0(snv1alt, snv2ref), paste0(snv1ref, snv2alt)),
                                 labels = c("Num_ref_ref", "Num_alt_alt", "Num_alt_ref", "Num_ref_alt"))))
    
    return(count_data)
  }
}

