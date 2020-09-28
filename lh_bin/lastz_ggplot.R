#!/usr/bin/env Rscript

library("ggplot2")

args = commandArgs(trailingOnly=TRUE)

pdf(paste(args[1], ".gg.pdf", sep=""))

fname = args[1]

test<-read.table(fname, header=TRUE,comment.char = "")

forward <- subset(test, test$strand2 == '+')
backward <- subset(test, test$strand2 == '-')


b <- ggplot(test, aes(x = zstart1, y = zstart2.))

b + geom_segment(aes(xend=end1, yend=end2., color = strand2), size=1)

##plot(NA, xlim=c(0,max(test$end1)), ylim=c(0,max(test$end2.)),
##     xlab="ref", ylab="qry", main=fname)
##
##segments(forward$zstart1,forward$zstart2.,forward$end1, forward$end2., col='blue', lwd=1) 
##segments(backward$zstart1,backward$zstart2.,backward$end1, backward$end2., col='red', lwd=1 )

dev.off()
