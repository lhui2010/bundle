#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output_file <- paste(args[1], "pdf", sep=".")
pdf(output_file, width=7, height=28)

library(ggplot2)

data=read.table(args[1], header=T)

library(ggplot2)
library(gggenes)

#setwd("/Users/liuhui/OneDrive/LH_KIB/Research/LegumePhylogenomics/PhyLegR/")
data$molecule <- factor(data$molecule, levels = unique(data$molecule))
#this_table$REF = factor(this_table$REF, levels = unique(this_table$REF))
ggplot(data, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow(arrowhead_width=grid::unit(0, "mm"),
                                    arrowhead_height = grid::unit(0, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
    scale_fill_brewer(palette = "Set3")  + theme_genes()
  dev.off()
