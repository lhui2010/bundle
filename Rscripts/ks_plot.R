#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

self_kaks = read.table(args[1], header = T)
pdf_name = paste(args[1], 'pdf', sep='.')

library(ggplot2)
self_kaks_filter = self_kaks[self_kaks$Ks<6,]

# a=data.frame(Ks=self_kaks_filter$Ks, Type=self_kaks_filter$file)

p <- ggplot(self_kaks_filter, aes(Ks, after_stat(density),  linetype="a")) +
    geom_density(bw = 0.05,alpha=0.5, position="identity")  #+
#p <- ggplot(self_kaks_filter, aes(Ks, after_stat(density), color = file, linetype="a")) +
#     facet_wrap(vars(Fam))
#    facet_grid(rows=vars(Fam))
ggsave(pdf_name, p, width=16, height=9)

