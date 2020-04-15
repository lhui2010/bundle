#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

pdf(paste(args[1], ".pdf", sep=""))

fname = args[1]

test<-read.table(fname, header=TRUE)

forward <- subset(test, test$strand2 == '+')
backward <- subset(test, test$strand2 == '-')

plot(NA, xlim=c(0,max(test$end1)), ylim=c(0,max(test$end2.)),
     xlab="ref", ylab="qry", main=fname)

segments(forward$zstart1,forward$zstart2.,forward$end1, forward$end2., col='blue') 
segments(backward$zstart1,backward$zstart2.,backward$end1, backward$end2., col='red' )

dev.off()
