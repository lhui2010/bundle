#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output_file <- paste("LTR_insert_fluctruation.pdf", sep="")
pdf(output_file, width=12, height=8)
par(mfrow=c(2,1))

title<-c("upstream", "downstream")

cutoff = 10000

for (i in 2:3){
    input_file <- args[1]
    sp_name <- title[i-1]
#output_file <- paste(input_file, "hist.pdf", sep="")
    a<-read.table(input_file)
    #hist(a[,3], breaks=seq(0,7.5,0.1), col='blue', main = sp_name,
#     ylab="Count of genes", ylim = c(0,8000), xaxt='n',xlab="LTR insert variations")
    b<-abs(a[,i])
    hist(b[b<cutoff & b>0], breaks=seq(0,cutoff,cutoff/100), col='blue', main = sp_name,
     ylab="Count of genes", xlab="LTR insert variations" )
#    axis(side =1, at=seq(0,7,1), labels=seq(0,7,1))
#    axis(side =1, at=seq(0,7,1), labels=seq(0,7,1))

}

dev.off()

