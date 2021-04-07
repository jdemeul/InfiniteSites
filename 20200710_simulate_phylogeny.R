### simulate sequences and tree

library(rtracklayer)
library(Biostrings)

cstring <- DNAStringSet(paste0(rep("C", 1000), collapse = ""))
cstring

ncycles <- 10

# 2^c(0:10)

dnastringlist <- list()
mutposlist <- list()
mutposlist[[1]] <- integer()
dnastringlist[[1]] <- cstring

for (cycleidx in 1:ncycles) {
# cycleidx <- 1
print(paste0("This is cycle ", cycleidx))

# take previous generation and create children

parents <- dnastringlist[[cycleidx]]
children <- parents[rep(1:length(parents), rep(2, length(parents)))]

mutantpositions <- sapply(X = vmatchPattern(pattern = "C", subject = children), FUN = function(x) sample(x = start(x), size = 1))
mutantpositions <- split(x = IRanges(start = mutantpositions, width = 1), f = 1:length(mutantpositions))
mutposlist[[cycleidx+1]] <- mutantpositions
dnastringlist[[cycleidx+1]] <- replaceAt(x = children, at = mutantpositions, value = "T")

}


mutposlist[[11]]
print(as.character(dnastringlist[[11]][1]))
print(as.character(dnastringlist[[11]][2]))


##########
# single genome

mutloads <- c(1e3, 2.5e3, 5e3, 1e4, 2.5e4, 5e4, 1e5)
genomesize <- 3e9/6

nlineages <- 2^c(0:6)
# mutlist <- list()
# plotlist <- list()
# 
nsims <- 100
# 
# for (mutload in mutloads) {
#   for (n in nlineages) {
#     for (sim in 1:nsims) {
#       mutlist[[paste0("lin", n)]][sim] <- sum(duplicated(1 + ( abs(sample(x = (-n*genomesize):(n*genomesize), size = n*mutload, replace = F)) %% genomesize )))
#     }
#   }
#   mutlistdf <- melt(data.frame(t(do.call(rbind, mutlist))))
#   mutlistdf$variable <- as.integer(gsub(pattern = "lin", replacement = "", x = mutlistdf$variable))
#   plotlist[paste0("mutload_", mutload)] <- mutlistdf
#   rm(mutlistdf)
# }

library(reshape2)

dosim <- function(mutload, nlineages, nsims, genomesize) {
  mutlist <- list()
  for (n in nlineages) {
    for (sim in 1:nsims) {
      mutlist[[paste0("lin", n)]][sim] <- sum(duplicated(1 + ( abs(sample(x = (-n*genomesize):(n*genomesize), size = n*mutload, replace = F)) %% genomesize )))
    }
  }
  mutlistdf <- melt(data.frame(t(do.call(rbind, mutlist))))
  mutlistdf$variable <- as.integer(gsub(pattern = "lin", replacement = "", x = mutlistdf$variable))
  return(mutlistdf)
}

# temp <- dosim(mutload = 1.5e4, nlineages = c(1,2,10), nsims = 10, genomesize = 3e9/6)
# temp

library(parallel)
outmuts <- mclapply(X = mutloads, FUN = dosim, nlineages = nlineages, nsims = nsims, genomesize = genomesize, mc.preschedule = T, mc.cores = 10)
outmutsbound <- data.frame(do.call(rbind, outmuts))
outmutsbound$mutload <- paste0("ml", rep(x = mutloads, times = rep(nsims*length(nlineages), times = length(mutloads)))/1e3)
outmutsbound$mutload <- factor(outmutsbound$mutload, levels = paste0("ml", mutloads/1e3))

library(ggplot2)


p1 <- ggplot(data = outmutsbound, mapping = aes(x = variable, y = value + 1, group = mutload, colour = mutload)) + 
  geom_point(shape = NA, position = position_jitter(width = 0.01)) +
  geom_smooth() + scale_y_log10() + scale_x_log10() + scale_colour_brewer(type = "seq", palette = "Oranges") +
  theme_minimal() + annotation_logticks(sides = "lb") + labs(x = "# samples", y = "# infinite sites violations + 1")
p1
ggsave(plot = p1, filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/figures/multisamplesimulator_violationcounts.pdf", width = 15, height = 12)
saveRDS(object = outmutsbound, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/multisamplesimulator_violationcounts.RDS")



####
mutload <- 1e5
genomesizes <- 3e9/c(2^(c(0:6)))

outmat <- t(replicate(n = 100,expr = sapply(X = genomesizes, FUN = function(x) sum(duplicated( abs(sample(x = -x:x, size = mutload, replace = F)))))))
colMeans(outmat)
plot(log10(genomesizes), log10(colMeans(outmat)))

(log10(genomesizes[length(genomesizes)]) - log10(genomesizes[1]))/(log10(colMeans(outmat)[length(genomesizes)]) - log10(colMeans(outmat)[1]))

nviol <- sum(duplicated(1 + ( abs(sample(x = (-n*genomesize):(n*genomesize), size = n*mutload, replace = F)) %% genomesize )))
