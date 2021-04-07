### comparison between two simulators for paper


ssize <- 1e4
checkparevo <- function(ssize) {
  v1 <- sample(x = -3e9:2e9, size = ssize, replace = F)
  v2 <- sample(x = -3e9:3e9, size = ssize, replace = F)
  v3 <- sample(x = -3e9:2e9, size = ssize, replace = F)
  v4 <- sample(x = -3e9:3e9, size = ssize, replace = F)
  nbial1 <- sum(duplicated(abs(v1)))
  nbial2 <- sum(duplicated(abs(v2)))
  nbial3 <- sum(duplicated(abs(v1)))
  nbial4 <- sum(duplicated(abs(v2)))
  npar2 <- sum(duplicated(abs(c(v1,v2))))
  npar3 <- sum(duplicated(abs(c(v1,v2,v3))))
  npar4 <- sum(duplicated(abs(c(v1,v2,v3,v4))))
  return(c(nbial1, nbial2, npar2, npar3, npar4))
}

library(parallel)
# outtemp <- lapply(2^(0:21), FUN = checkparevo)
outtemp <- mclapply(X = c(seq(from = 1e4, to = 1e5, by = 1e4), seq(from = 1e5, to = 1e6, by = 2.5e4), seq(from = 1e6, to = 2e6, by = 1e5)),
                    FUN = checkparevo, mc.cores = 10)
outtemp <- data.frame(do.call(rbind, outtemp))
outtemp$n <- c(seq(from = 1e4, to = 1e5, by = 1e4), seq(from = 1e5, to = 1e6, by = 2.5e4), seq(from = 1e6, to = 2e6, by = 1e5))
outtemp
p1 <- ggplot(data = outtemp, mapping = aes(x = n)) + geom_point(mapping = aes(y = X1), color = "blue") + geom_point(mapping = aes(y = X3), color = "red") + geom_point(mapping = aes(y = X4), color = "pink") + geom_point(mapping = aes(y = X5), color = "green") + scale_y_log10() + scale_x_log10() + annotation_logticks()
p1

colMeans((outtemp/(outtemp$X1*4))[-(1:20),])
