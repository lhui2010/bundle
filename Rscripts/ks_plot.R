#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

self_kaks = read.table(args[1], header = T)
pdf_name = paste(args[1], 'pdf', sep='.')

library(ggplot2)
self_kaks_filter = self_kaks[self_kaks$Ks<6,]

a=data.frame(Ks=self_kaks_filter$Ks, Type=self_kaks_filter$file)

p <- ggplot(a, aes(Ks, after_stat(density), color = Type, linetype="a")) +
    geom_density(bw = 0.05,alpha=0.5, position="identity")
ggsave(pdf_name, p)

