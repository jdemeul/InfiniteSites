#### clustering by mutation spectrum

bialsummary <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesBiallelicM2recall_summary.txt", as.is = T)
alsummary <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/InfSitesByAF_alphapt01_hetonly_summary_cleanhits.txt", as.is = T)

mutsigs <- read.delim(file = "/srv/shared/vanloo/ICGC_signatures/20180322_release/PCAWG_sigProfiler_SBS_signatures_in_samples_waliqID.csv", sep = ",", as.is = T)
mutsigs[, 4:68] <- mutsigs[, 4:68]/rowSums(mutsigs[, 4:68])

mutsigssub <- mutsigs[mutsigs$Sample.Name %in% c(bialsummary[bialsummary$nbiallelics >= 1, "sampleid"], alsummary[alsummary$nparallel >= 1,  "sampleid"]), ]
mutsigssub <- mutsigssub[ , c("Cancer.Types", "Sample.Name", "Accuracy", names(which(colSums(mutsigssub[, 4:68]) > 0)))]
mutsigpca <- prcomp(mutsigssub[, 4:ncol(mutsigssub)], center = T)
plot(mutsigpca)

mutsigssub2 <- cbind(mutsigssub[, 1:3], mutsigpca$x)

library(ggplot2)

p1 <- ggplot(data = mutsigssub2, mapping = aes(x = PC2, y = PC3)) + geom_point(mapping = aes(colour = Cancer.Types))
p1

temp <- kmeans(x = mutsigssub[, 4:ncol(mutsigssub)], centers = 4)
plot(temp)

library(Rtsne)
temp <- Rtsne(X = mutsigssub[, 4:ncol(mutsigssub)])
mutsigssub2 <- cbind(mutsigssub[, 1:3], temp$Y)
colnames(mutsigssub2)[4:5] <- c("tsne1", "tsne2")

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(mutsigssub2$Cancer.Types))), scheme = "tumour.subtype")
names(cvect) <- unique(mutsigssub2$Cancer.Types)
# cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", '#FF4500', '#FF4500')


p1 <- ggplot(data = mutsigssub2, mapping = aes(x = tsne1, y = tsne2)) + geom_point(mapping = aes(colour = Cancer.Types))
p1 <- p1 + scale_color_manual(values = cvect)
p1



#### use sims
simtots <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/results/isabreakdown_1plus1_only_writefrac/", pattern = "_infsites_permut_effgenfrac0.1.txt", recursive = T, full.names = T)
simtotsall <- lapply(simtots, FUN = read.delim, as.is = T)
simtotsall <- as.data.frame(do.call(rbind, lapply(simtotsall, function(x) x$counts_tot)))
simtotsall <- simtotsall/rowSums(simtotsall)
simtotsall$sampleid <- gsub(pattern = "_infsites_permut_effgenfrac0.1.txt", replacement = "", x = basename(simtots))
simtotsall$Cancer.Types <- mutsigs[match(x = simtotsall$sampleid, table = mutsigs$Sample.Name), "Cancer.Types"]
simtotsall <- simtotsall[!is.nan(rowSums(simtotsall[, 1:456])), ]

temp <- Rtsne(X = simtotsall[, 1:456])
simtotsall2 <- cbind(simtotsall[, 457:458], temp$Y)
colnames(simtotsall2)[3:4] <- c("tsne1", "tsne2")

cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(simtotsall2$Cancer.Types))), scheme = "tumour.subtype")
names(cvect) <- unique(simtotsall2$Cancer.Types)
# cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", '#FF4500', '#FF4500')


p1 <- ggplot(data = simtotsall2, mapping = aes(x = tsne1, y = tsne2)) + geom_point(mapping = aes(colour = Cancer.Types))
p1 <- p1 + scale_color_manual(values = cvect)
p1
