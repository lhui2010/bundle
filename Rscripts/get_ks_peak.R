#!/usr/bin/env Rscript

require(mclust)
#get_peak.R kaks_file peak
args = commandArgs(trailingOnly=TRUE)
self_kaks = read.table(args[1], header = T)
self_kaks_filter = self_kaks[self_kaks$Ks<6,]
X <- na.omit(self_kaks_filter$Ks)
# BIC <- mclustBIC(X)
# mod1 <- Mclust(X, x = BIC)
# summary(mod1, parameters = TRUE)
if(length(args) == 2)
{
    peak = args[2]
    mod4 <- densityMclust(X, G=peak)
    pdf(paste(args[1], '.kspeak.G', peak, '.pdf', sep=''))
} else
{
    mod4 <- densityMclust(X)
    pdf(paste(args[1], 'kspeak.Gnull', 'pdf', sep='.'))
}
print(mod4$parameters)
print(summary(mod4))
#plot(mod4, what="BIC")
plot(mod4, what = "density", data = X, breaks = 50)
