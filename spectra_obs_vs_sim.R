### linking observed spectra to predicted spectra (1+1 regions) both parallel & third alleles
library(ggplot2)
library(VariantAnnotation)
library(BSgenome.Hsapiens.1000genomes.hs37d5)




#### fions

generate_bases_types_trinuc <- function() {
  # generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  types <- c("A", "G", "T", "A", "C", "G")
  types2 <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(trinucleotides, ">", rep(types, rep(16,6))))
  return(list(bases = bases, types = types, types2 = types2, trinucleotides = trinucleotides,
              trinucleotides_mutations = trinucleotides_mutations))
}


get_mutspectrum <- function(mutations, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5) {
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  if (!"trinuc" %in% colnames(mutations))
    mutations$trinuc <- as.character(get_trinuc_context(mutations = mutations, bsgenome = bsgenome))
  
  # reverse complement and data augmentation
  revcomp <- data.frame(alt = as.character(complement(DNAStringSet(mutations$alt))),
                        trinuc = as.character(reverseComplement(DNAStringSet(mutations$trinuc))),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                            paste0(mutations$trinuc, ">", mutations$alt),
                            paste0(revcomp$trinuc, ">", revcomp$alt)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations"]])
  
  # compute frequencies and normalised probabilities
  muttype_freq <- c(table(mut_full))
  return(muttype_freq)
}



get_mutspectrum_biallelic <- function(mutations, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5) {
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  if (!"trinuc" %in% colnames(mutations))
    mutations$trinuc <- as.character(get_trinuc_context(mutations = mutations, bsgenome = bsgenome))
  
  altbases <- split(cbind(mutations$alt1, mutations$alt2), f = 1:nrow(mutations))
  alt <- ifelse(mutations$ref == "C", ifelse(sapply(altbases, setequal, y = c("A","G")), "A+G", 
                                             ifelse(sapply(altbases, setequal, y = c("T","G")), "G+T", "T+A")),
                ifelse(mutations$ref == "G", ifelse(sapply(altbases, setequal, y = c("T","C")), "A+G", 
                                                    ifelse(sapply(altbases, setequal, y = c("A","C")), "G+T", "T+A")),
                       ifelse(mutations$ref == "T", ifelse(sapply(altbases, setequal, y = c("A","C")), "A+C", 
                                                           ifelse(sapply(altbases, setequal, y =c("C","G")), "C+G", "G+A")),
                              ifelse(sapply(altbases, setequal, y = c("T","G")), "A+C", 
                                     ifelse(sapply(altbases, setequal, y =c("G","C")), "C+G", "G+A")))
                ))
  
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"),
                            paste0(mutations$trinuc, ">", alt), 
                            paste0(as.character(reverseComplement(DNAStringSet(mutations$trinuc))), ">", alt)),
                     levels = paste0(base_type_trinuc_info[["trinucleotides_mutations"]], "+", rep(c("G", "T", "A", "C", "G", "A"), rep(16,6))))
  
  # compute frequencies and normalised probabilities
  muttype_freq <- c(table(mut_full))
  return(muttype_freq)
}



# get the trinucleotide context of a set of point mutations
get_trinuc_context <- function(mutations, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, size = 1) {
  sequences <- getSeq(x = bsgenome, names = mutations$chr, start = mutations$start - size, end = mutations$end + size, strand = "+")
  return(sequences)
}



# get the trinucleotide content of the genome
get_trinuc_normalisation_factors <- function(bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5) {
  reference_trinucleotides <- unique(generate_bases_types_trinuc()$trinucleotides)
  nonreference_trinucleotides <- as.character(reverseComplement(DNAStringSet(reference_trinucleotides)))
  
  trinuc_counts_all <- colSums(trinucleotideFrequency(x = getSeq(x = bsgenome), as.prob = F))
  trinuc_counts <- trinuc_counts_all[reference_trinucleotides] + trinuc_counts_all[nonreference_trinucleotides]
  trinuc_counts / sum(trinuc_counts)
  
  return(trinuc_counts)
}



### compare spectra by taking simulated as exact and deriving HDI of observed using Posterior Dirichlet Multinomial
annotate_spectrum <- function(sampleid, trinucs, simsubset = "1plus1") {
  
  print(paste0("Running sample ", sampleid))
  
  bases <- generate_bases_types_trinuc()
  
  # load all SNVs
  snvfile <- list.files(path = "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/", pattern = paste0(sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz$"), recursive = T, full.names = T)
  # snvfile <- paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/snv_mnv/", sampleid, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
  snvs <- readVcf(file = snvfile, genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
  snvs <- data.frame(chr = seqnames(snvs), start = start(snvs), end = end(snvs), ref = as.character(ref(snvs)), alt = as.character(unlist(alt(snvs))))
  
  # load all parallel mutations (i.e. in all but LOH regions)
  # vafhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_snv_mnv_infSites_finalhits.txt"))
  vafhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_snv_mnv_infSites_finalhits.txt"))
  vafhitsdf <- read.delim(file = vafhitsfile, as.is = T)
  vafhitsdf$chr <- as.character(vafhitsdf$chr)
  # vafhitsdf$pos <- vafhitsdf$start
  
  ### read total simulated in het regions
  simfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_", simsubset, "_only_writefrac/", sampleid, "/", paste0(sampleid, "_infsites_permut_effgenfrac0.1.txt"))
  if (!file.exists(simfile)) {
    return(NULL)
  }
  simdf <- read.delim(file = simfile, as.is = T)
  # simdf$prob <- (1+simdf$counts_tot)/sum(simdf$counts_tot+1)
  
  # get full mutspectrum (background)
  allmutstab <- get_mutspectrum(mutations = snvs, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)

  plotdf <- simdf[simdf$type == "parallel", ]
  if (nrow(vafhitsdf) > 0 & sum(plotdf$counts_tot) > 0) {
    vafhitsdf$ref <- gsub(pattern = "RUE", replacement = "", x = as.character(vafhitsdf$ref))
    vafhitsdf$alt <- gsub(pattern = "RUE", replacement = "", x = as.character(vafhitsdf$alt))
    # bases <- generate_bases_types_trinuc()
    parmutstab <- get_mutspectrum(mutations = vafhitsdf, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
    
    plotdf$allobs <- allmutstab[plotdf$mut1]
    plotdf$parobs <- parmutstab[plotdf$mut1]
    
    plotdf$trinuc <- substr(plotdf$mut1, start = 1, stop = 3)
    
    # outprobs <- t(apply(X = MCMCpack::rdirichlet(n = 10000, alpha = plotdf$parobs+plotdf$allobs*96/sum(plotdf$allobs)),
    #                    MARGIN = 2, FUN = quantile, probs = c(0.025, .5, 0.975)))
    # outprobs <- as.data.frame(outprobs)
    # colnames(outprobs) <- c("pobs_low", "pobs_med", "pobs_hi")
    rmat <- MCMCpack::rdirichlet(n = 10000, alpha = plotdf$parobs+9.6*plotdf$allobs/sum(plotdf$allobs))
    # rmat <- MCMCpack::rdirichlet(n = 10000, alpha = plotdf$parobs+1)
    # rmat <- MCMCpack::rdirichlet(n = 10000, alpha = plotdf$parobs+1/96)
    # pmat <- t(apply(X = pmat, MARGIN = 1, FUN = rmultinom, n = 1, size = sum(plotdf$parobs)))
    rmat <- sweep(x = rmat, MARGIN = 2, trinucs[plotdf$trinuc], FUN = "/")
    rmat <- sweep(x = rmat, MARGIN = 1, rowSums(rmat), FUN = "/")
    
    outprobs <- apply(X = rmat, MARGIN = 2, FUN = quantile, probs = c(.025, .5, .975))
    outprobs <- as.data.frame(t(outprobs))
    colnames(outprobs) <- c("pobs_low", "pobs_med", "pobs_hi")
    
    plotdf <- cbind(plotdf, outprobs)
    
    
    plotdf$allobs_norm <- plotdf$allobs/trinucs[plotdf$trinuc]
    plotdf$allobs_norm <- plotdf$allobs_norm/sum(plotdf$allobs_norm)
    # plotdf$allobs_norm <- plotdf$allobs/sum(plotdf$allobs)
    
    plotdf$sim_norm <- plotdf$counts_tot/trinucs[plotdf$trinuc]
    plotdf$sim_norm <- plotdf$sim_norm/sum(plotdf$sim_norm)
    # plotdf$sim_norm <- plotdf$counts_tot/sum(plotdf$counts_tot)
    
    #### checking log odds or cosine similarity between the vectors
    
    # logpsim <- dmultinom(x = plotdf$parobs, prob = (plotdf$counts_tot+1/96)/sum(plotdf$counts_tot+1/96), log = T)
    # logpall <- dmultinom(x = plotdf$parobs, prob = (plotdf$allobs+1/96)/sum(plotdf$allobs+1/96), log = T)
    # cossim <- sum(plotdf$parobs*plotdf$counts_tot)/(sqrt(sum(plotdf$parobs^2))*sqrt(sum(plotdf$counts_tot^2)))
    # cosall <- sum(as.numeric(plotdf$parobs)*as.numeric(plotdf$allobs))/(sqrt(sum(plotdf$parobs^2))*sqrt(sum(plotdf$allobs^2)))
    
    parcossim <- quantile(apply(X = rmat, MARGIN = 1, FUN = function(x, y) sum(x*y)/(sqrt(sum(x^2))*sqrt(sum(y^2))), y = plotdf$sim_norm), probs = c(.025, .5, .975), na.rm = T)
    
    # cosall <- sum(sum(as.numeric(plotdf2$thirdobs)*as.numeric(plotdf2$allobs)))/(sqrt(sum(plotdf2$thirdobs^2))*sqrt(sum(plotdf2$allobs^2)))
    parcosall <- quantile(apply(X = rmat, MARGIN = 1, FUN = function(x, y) sum(x*y)/(sqrt(sum(x^2))*sqrt(sum(y^2))), y = plotdf$allobs_norm), probs = c(.025, .5, .975), na.rm = T)
    
    rm(rmat)
    
    
    # apply(X = pmat, MARGIN = 1, FUN = function(x) x)
    
    plotdf$type <- paste0(substr(plotdf$mut1,2,2), ">", substr(plotdf$mut1,5,5))
    plotdf$mut1 <- factor(plotdf$mut1, levels = bases$trinucleotides_mutations)
    
    # if (T) {
    #   p1 <- ggplot(data = plotdf, mapping = aes(x = mut1)) + geom_col(mapping = aes(y = allobs_norm, fill = type), show.legend = F, alpha = .6)
    #   p1 <- p1 + geom_segment(data = data.frame(start = seq(0.65, 96, 16), end = seq(16.35, 97, 16), type = bases$types2),
    #                           mapping = aes(x = start, xend = end, y = 1.25*max(plotdf$allobs_norm), yend = 1.25*max(plotdf$allobs_norm), color = type), size = 3, show.legend = F, alpha = .6)
    #   p1 <- p1 + geom_text(data = data.frame(type = bases$types2, pos = seq(8, 96, 16)),
    #                        mapping = aes(x = pos, y = 1.3*max(plotdf$allobs_norm), color = type, label = type),
    #                        show.legend = F, family = "mono")
    #   p1 <- p1 + theme_bw() + scale_x_discrete(labels = plotdf$trinuc) + scale_fill_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
    #     scale_color_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
    #     coord_cartesian(ylim = c(0, 1.15*max(plotdf$allobs_norm)), clip="off") + ylab(label = "Probability") +
    #     theme(axis.text.x = element_text(angle = 90, family = "mono"), panel.grid.major.x = element_blank(),
    #           axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(2,.5,.5,.5), "lines"))
    #   ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", sampleid, "_spectrum.pdf"), plot = p1, width = 16, height = 4)
    # }
    
    p1 <- ggplot(data = plotdf, mapping = aes(x = mut1)) + 
      geom_col(mapping = aes(y = allobs_norm, fill = type), show.legend = F, alpha = .6) +
      geom_pointrange(mapping = aes(y = pobs_med, ymin = pobs_low, ymax = pobs_hi)) +
      geom_point(mapping = aes(y = sim_norm), colour = "#ff7f00", alpha = .8)
    p1 <- p1 + geom_segment(data = data.frame(start = seq(0.65, 96, 16), end = seq(16.35, 97, 16), type = bases$types2),
                            mapping = aes(x = start, xend = end, y = 1.25*max(plotdf$pobs_hi), yend = 1.25*max(plotdf$pobs_hi), color = type), size = 3, show.legend = F, alpha = .6)
    p1 <- p1 + geom_text(data = data.frame(type = bases$types2, pos = seq(8, 96, 16)),
                         mapping = aes(x = pos, y = 1.3*max(plotdf$pobs_hi), color = type, label = type),
                         show.legend = F, family = "mono")
    p1 <- p1 + theme_bw() + scale_x_discrete(labels = plotdf$trinuc) + scale_fill_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
      scale_color_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
      coord_cartesian(ylim = c(0, 1.15*max(plotdf$pobs_hi)), clip="off") + ylab(label = "Probability") +
      theme(axis.text.x = element_text(angle = 90, family = "mono"), panel.grid.major.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(2,.5,.5,.5), "lines"))
    # p1
    
    outfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_parallel_snv_spectra_", simsubset))
    # outfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_parallel_snv_spectra_", simsubset))
    ggsave(filename = paste0(outfile, ".pdf"), plot = p1, width = 16, height = 4, useDingbats=FALSE)
    write.table(x = plotdf, file = paste0(outfile, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  } else {
    parcosall <- as.numeric(rep(NA, 3))
    parcossim <- as.numeric(rep(NA, 3))
  }


  
  ### third allele stuff
  # thirdhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants.txt")
  thirdhitsfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_variants_20210217.txt")
  thirdhitsdf <- read.delim(file = thirdhitsfile, as.is = T)

  # if (simsubset == "1plus1") {
  #   thirdhitsdf <- thirdhitsdf[which(thirdhitsdf$major_cn == 1 & thirdhitsdf$minor_cn == 1), ]
  # }
  
  thirdhitsdf <- thirdhitsdf[thirdhitsdf$sampleid == sampleid, ]
  colnames(thirdhitsdf)[c(1, 7)] <- c("chr", "ref")
  thirdhitsdf$ref <- gsub(pattern = "RUE", replacement = "", x = as.character(thirdhitsdf$ref))
  thirdhitsdf$alt1 <- gsub(pattern = "RUE", replacement = "", x = as.character(thirdhitsdf$alt1))
  thirdhitsdf$alt2 <- gsub(pattern = "RUE", replacement = "", x = as.character(thirdhitsdf$alt2))
  
  plotdf2 <- simdf[simdf$type == "third_allele", ]
  if (nrow(thirdhitsdf) > 0 & sum(plotdf2$counts_tot) > 0) {
      
    b1 <- split(cbind(substr(plotdf2$mut1,5,5), substr(plotdf2$mut2,5,5)), f = 1:nrow(plotdf2))
    b2 <- substr(plotdf2$mut1,2,2)
    plotdf2$type <- paste0(b2, ">", ifelse(b2 == "C", ifelse(sapply(b1, setequal, y = c("A","G")), "A+G", 
                                                             ifelse(sapply(b1, setequal, y = c("T","G")), "G+T", "T+A")),
                                           ifelse(sapply(b1, setequal, y = c("A","C")), "A+C", 
                                                  ifelse(sapply(b1, setequal, y =c("C","G")), "C+G", "G+A"))))
    plotdf2$fulltype <- paste0(substr(plotdf2$mut1, start = 1, stop = 3), substr(plotdf2$type, start = 2, stop = 5))
    
    # plotdf2$primary <- paste0(substr(plotdf2$mut1, start = 1, stop = 4), substr(plotdf2$type, start = 3, stop = 3))
    plotdf2$trinuc <- substr(plotdf2$mut1, start = 1, stop = 3)
    
    # thirdsimdf$mut1 <- factor(thirdsimdf$mut1, levels = bases$trinucleotides_mutations)
    # vafhitsdf$pos <- vafhitsdf$start
    
    plotdf2$allobs <- allmutstab[substr(plotdf2$fulltype, 1,5)]
    # plotdf2$allobsOR <- plotdf2$allobs + allmutstab[paste0(substr(plotdf2$fulltype, 1,4), substr(plotdf2$fulltype, 7,7))]
    
    # debug(get_mutspectrum_biallelic)
    plotdf2$fulltype <- factor(x = plotdf2$fulltype, levels = paste0(bases[["trinucleotides_mutations"]], "+", rep(c("G", "T", "A", "C", "G", "A"), rep(16,6))))
    plotdf2 <- plotdf2[order(plotdf2$fulltype, decreasing = F), ]
    thirdmutstab <- get_mutspectrum_biallelic(mutations = thirdhitsdf, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
    plotdf2$thirdobs <- thirdmutstab[as.character(plotdf2$fulltype)]
    
    # temp <- boot(data = thirdhitsdf, statistic = function(df, i) {get_mutspectrum_biallelic(df[i,])}, R = 100)
    # pmat <- sweep(x = temp$t, MARGIN = 2, trinucs[plotdf2$trinuc], FUN = "/")
    # pmat <- sweep(x = temp$t, MARGIN = 1, rowSums(temp$t), FUN = "/")
    
    rmat <- MCMCpack::rdirichlet(n = 10000, alpha = plotdf2$thirdobs+1/96)
    # pmat <- t(apply(X = pmat, MARGIN = 1, FUN = rmultinom, n = 1, size = sum(plotdf2$thirdobs)))
    rmat <- sweep(x = rmat, MARGIN = 2, trinucs[plotdf2$trinuc], FUN = "/")
    rmat <- sweep(x = rmat, MARGIN = 1, rowSums(rmat), FUN = "/")
    
    outprobs <- apply(X = rmat, MARGIN = 2, FUN = quantile, probs = c(.025, .5, .975))
    outprobs <- as.data.frame(t(outprobs))
    colnames(outprobs) <- c("pobs_low", "pobs_med", "pobs_hi")
    
    plotdf2 <- cbind(plotdf2, outprobs)
    
    
    plotdf2$allobs_norm <- plotdf2$allobs/trinucs[plotdf2$trinuc]
    plotdf2$allobs_norm <- plotdf2$allobs_norm/sum(plotdf2$allobs_norm)
    # plotdf$allobs_norm <- plotdf$allobs/sum(plotdf$allobs)
    
    plotdf2$allobsOR_norm <- plotdf2$allobs_norm + plotdf2[match(table = substr(plotdf2$fulltype, 1, 5), x = paste0(substr(plotdf2$fulltype, 1,4), substr(plotdf2$fulltype, 7,7))), "allobs_norm"]
    # plotdf2$allobsOR_norm <- exp(log(plotdf2$allobsOR) - log(trinucs[plotdf2$trinuc]))
    # plotdf2$allobsOR_norm <- plotdf2$allobsOR_norm/sum(plotdf2$allobsOR_norm)
    
    plotdf2$sim_norm <- plotdf2$counts_tot/trinucs[plotdf2$trinuc]
    plotdf2$sim_norm <- plotdf2$sim_norm/sum(plotdf2$sim_norm)
    # plotdf$sim_norm <- plotdf$counts_tot/sum(plotdf$counts_tot)
    
    #### checking log odds or cosine similarity between the vectors
    plotdf2$typesplit <- substr(plotdf2$type, 1,3)
    plotdf3 <- plotdf2
    plotdf3$typesplit <- paste0(substr(plotdf3$type, 1,2), substr(plotdf3$type, 5,5))
    plotdf3$allobs_norm <- plotdf3[match(table = substr(plotdf3$fulltype, 1, 5), x = paste0(substr(plotdf3$fulltype, 1,4), substr(plotdf3$fulltype, 7,7))), "allobs_norm"]
    
    plotdf3 <- rbind(plotdf2, plotdf3)
    plotdf3$typesplit <- factor(plotdf3$typesplit, levels = c("C>A", "C>G", "C>T","T>A", "T>C","T>G"))
    
    # logpsim <- dmultinom(x = plotdf2$thirdobs, prob = (plotdf2$counts_tot+1/96)/sum(plotdf2$counts_tot+1/96), log = T)
    # logpall <- dmultinom(x = plotdf2$thirdobs, prob = (plotdf2$allobs+1/96)/sum(plotdf2$allobs+1/96), log = T)
    # cossim <- sum(plotdf2$thirdobs*plotdf2$counts_tot)/(sqrt(sum(plotdf2$thirdobs^2))*sqrt(sum(plotdf2$counts_tot^2)))
    thirdcossim <- quantile(apply(X = rmat, MARGIN = 1, FUN = function(x, y) sum(x*y)/(sqrt(sum(x^2))*sqrt(sum(y^2))), y = plotdf2$sim_norm), probs = c(.025, .5, .975), na.rm = T)
    
    # cosall <- sum(sum(as.numeric(plotdf2$thirdobs)*as.numeric(plotdf2$allobs)))/(sqrt(sum(plotdf2$thirdobs^2))*sqrt(sum(plotdf2$allobs^2)))
    thirdcosall <- quantile(apply(X = rmat, MARGIN = 1, FUN = function(x, y) sum(x*y)/(sqrt(sum(x^2))*sqrt(sum(y^2))), y = plotdf2$allobsOR_norm), probs = c(.025, .5, .975), na.rm = T)
    
    rm(rmat)
    
    p1 <- ggplot(data = plotdf2, mapping = aes(x = fulltype)) + 
      geom_col(data = plotdf3, mapping = aes(y = allobs_norm, fill = typesplit), position = "stack", show.legend = F, alpha = .6) +
      geom_pointrange(mapping = aes(y = pobs_med, ymin = pobs_low, ymax = pobs_hi)) + 
      geom_point(mapping = aes(y = sim_norm), colour = "#ff7f00", alpha = .8)
    p1 <- p1 + geom_segment(data = data.frame(start = seq(0.65, 96, 16), end = seq(16.35, 97, 16), type = paste0(bases$types2, "+", c("G", "T", "A", "C", "G", "A"))),
                            mapping = aes(x = start, xend = end, y = 1.25*max(plotdf2$pobs_hi), yend = 1.25*max(plotdf2$pobs_hi), color = type), size = 3, show.legend = F, alpha = .6)
    p1 <- p1 + geom_text(data = data.frame(type = paste0(bases$types2, "+", c("G", "T", "A", "C", "G", "A")), pos = seq(8, 96, 16)),
                         mapping = aes(x = pos, y = 1.3*max(plotdf2$pobs_hi), color = type, label = type),
                         show.legend = F, family = "mono")
    p1 <- p1 + theme_bw() + scale_x_discrete(labels = plotdf2$trinuc) + scale_fill_manual(values = c("C>A" = "#15bcee","C>G" = "#000000","C>T" = "#e32926","T>A" = "#999999","T>C" = "#a1ce63","T>G" = "#ebc6c4")) +
      scale_color_manual(values = c("C>A+G" = "#15bcee","C>G+T" = "#000000","C>T+A" = "#e32926","T>A+C" = "#999999","T>C+G" = "#a1ce63","T>G+A" = "#ebc6c4")) +
      coord_cartesian(ylim = c(0, 1.15*max(plotdf2$pobs_hi)), clip="off") + ylab(label = "Probability") +
      theme(axis.text.x = element_text(angle = 90, family = "mono"), panel.grid.major.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(2,.5,.5,.5), "lines"))
    # p1
    outfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_thirdallele_snv_spectra_", simsubset))
    # outfile <- file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_thirdallele_snv_spectra_", simsubset))
    ggsave(filename = paste0(outfile, ".pdf"), plot = p1, width = 16, height = 4, useDingbats=FALSE)
    write.table(x = plotdf2, file = paste0(outfile, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  } else {
    thirdcosall <- as.numeric(rep(NA, 3))
    thirdcossim <- as.numeric(rep(NA, 3))
  }
  
  cosdf <- as.data.frame(rbind(parcosall, parcossim, thirdcosall, thirdcossim))
  colnames(cosdf) <- c("cosine_low", "cosine_med", "cosine_hi")
  cosdf$type <- c("parallel","parallel","third_allele","third_allele")
  cosdf$comparison <- c("all","simulated","all","simulated")
  # write.table(x = cosdf, file = file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20190416_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_parallel_thirdallele_cosinesims_", simsubset, ".txt")), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = cosdf, file = file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/20210212_vafpipeline_out_alphapt1_hetonly/", sampleid, paste0(sampleid, "_parallel_thirdallele_cosinesims_", simsubset, ".txt")), quote = F, sep = "\t", row.names = F, col.names = T)
  
  return(NULL)
}

#### fions








##### checking mtuation spectra

# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/mutation_spectrum_analysis_functions_20181206.R")


# temp <- get_wg_trinuc_normalisation_factors(bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
# plot_mutationspectrum(mutations = snvsdf, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "all2", normalise = T)


# plot_mutationspectrum(mutations = snvsdf, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "squared", normalise = F)


# snvsdfsub <- vafhitsdf[, c("chr", "start", "ref", "alt")]
# colnames(snvsdfsub)[2] <- "pos"
# plot_mutationspectrum(mutations = snvsdfsub, trinuc_freq = temp, sample = sampleid, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5, aspdf = T, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/", suffix = "ISAviol2", normalise = T)


trinucs <- get_trinuc_normalisation_factors(bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
trinucs <- trinucs/sum(trinucs)


# sampleid <- "2df02f2b-9f1c-4249-b3b4-b03079cd97d9"
# sample id <- "deb9fbb6-656b-41ce-8299-554efc2379bd"
# sampleid <- "c9f91ded-3b04-4cd1-8ea6-bbc635a8a4f0"
# sampleid <- "b07bad52-d44c-4b27-900a-960985bfadec"
# sampleid <- "c9f91ded-3b04-4cd1-8ea6-bbc635a8a4f0"
# sampleid <- "760881cc-c623-11e3-bf01-24c6515278c0"
# sampleid <- "bd3e88b3-b37c-4641-85fa-d8125ba324ca"
# sampleid <- "3b20d548-2a7d-4031-85a1-425ca7201d7a"

bialsummary <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary_20210217.txt", as.is = T)
alsummary <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits_v2_2021.txt", as.is = T)

# rslurmdf <- data.frame(sampleid = union(bialsummary$sampleid[bialsummary$nbiallelics > 0], alsummary$sampleid), stringsAsFactors = F)
rslurmdf <- data.frame(sampleid = alsummary[which(alsummary$ndivergent >= 10 | alsummary$npar_vaf >= 10), "sampleid"], stringsAsFactors = F)

# debug(annotate_spectrum)
# annotate_spectrum(sampleid = sampleid, trinucs = trinucs, simsubset = "het")
# lapply(X = sampleid, FUN = annotate_spectrum, trinucs = trinucs, simsubset = "1plus1")
# mclapply(X = rslurmdf$sampleid, FUN = annotate_spectrum, trinucs = trinucs, simsubset = "het", mc.preschedule = T, mc.cores = 8)
lapply(X = rslurmdf$sampleid, FUN = annotate_spectrum, trinucs = trinucs, simsubset = "het")


