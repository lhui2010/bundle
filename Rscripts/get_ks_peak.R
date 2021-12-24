#!/usr/bin/env Rscript

require(mclust)
args = commandArgs(trailingOnly=TRUE)
self_kaks = read.table(args[1], header = T)
self_kaks_filter = self_kaks[self_kaks$Ks<6,]
X <- na.omit(self_kaks_filter$Ks)
# BIC <- mclustBIC(X)
# mod1 <- Mclust(X, x = BIC)
# summary(mod1, parameters = TRUE)
mod4 <- densityMclust(X)
print(mod4$parameters)
print(summary(mod4))
pdf(paste(args[1], 'kspeak.pdf', sep='.'))
#plot(mod4, what="BIC")
plot(mod4, what = "density", data = X, breaks = 50)
