#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output_file <- paste(args[1], "pdf", sep=".")
pdf(output_file, width=12, height=8)
par(mfrow=c(2,1))

input_file <- args[1]
sp_name <- args[1]
a<-read.table(input_file, header = T)
b<-a[,1]
cutoff = max(b)
min = 0.5
hist(b[b<=cutoff & b>=min], breaks=seq(min,cutoff,(cutoff-min)/100), col='blue', main = sp_name,
 ylab="Count", xlab="Observed value" )

boxplot(b[b<=cutoff & b>=min])


dev.off()

