#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output_file <- paste(args[1], "pdf", sep=".")
pdf(output_file, width=12, height=8)

library(ggplot2)

data=read.table(args[1], header=T)

#head(data)

max_ks=6

data <- data$Ks[data$Ks <= max_ks];
data <- data[data >= 0.1];
data <- na.omit(data)

foo=args[1]

#     main=expression(paste('K'[s], ' Plot for', foo)),
hist(data, breaks=seq(0, max_ks, by=0.05),
     main=paste('Ks', ' Plot for', foo),
     xlab=expression(paste('Pairwise', ' K'[s])), axes=T);

# ggplot(b, aes(V1)) + geom_histogram(binwidth=0.1)


