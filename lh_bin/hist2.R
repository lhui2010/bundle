#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output_file <- paste(args[1], "pdf", sep=".")
pdf(output_file, width=12, height=8)
par(mfrow=c(2,1))

cutoff = 38000

input_file <- args[1]
sp_name <- args[1]
a<-read.table(input_file)
b<-abs(a[,3] - a[,2])
hist(b[b<cutoff & b>0], breaks=seq(0,cutoff,cutoff/100), col='blue', main = sp_name,
 ylab="Count of genes", xlab="LTR insert variations" )

boxplot(b[b<cutoff])


dev.off()

