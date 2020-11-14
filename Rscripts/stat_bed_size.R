#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
a=read.table(args[1], header = F)
#out = data.frame(rep("", 10))
cat("ChromosomeID\tMean\tMax\tNumber\tSize\n")
#for (i in 1:10){
for (i in unique(a$V1)){
    a_subset = subset(a, a$V1==i)
    len_list = a_subset$V3 - a_subset$V2
    cat(paste(i, "\t", as.integer(mean(len_list)), 
    max(len_list), 
    length(len_list), 
    sum(len_list),
    "\n")) 
#    print(paste(i, "\t", mean(len_list), max(len_list), length(len_list)), 
#    quote = F
#    )
#    print(summary(a_subset$V3 - a_subset$V2))
#    print(paste("Number: ", length(a_subset$V3 - a_subset$V2)))
}
cat(paste("Total\t", 
as.integer(mean(a$V3 - a$V2)), 
max(a$V3-a$V2), 
length(a$V2), 
sum(a$V3-a$V2),
"\n"))
