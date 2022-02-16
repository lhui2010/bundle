#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output_file <- paste(args[1], "pdf", sep=".")
pdf(output_file, width=12, height=8)

library(ggplot2)

a=read.table(args[1], header=F)

ggplot(a, aes(V1)) + geom_histogram()


