#!/usr/bin/env Rscript
#theme_Publication <- function(base_size=14, base_family="helvetica") {
theme_Publication <- function(base_size=14, base_family="arial") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
library(ggplot2)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)

#tb <- read.table("chr02.add_ref.OsSPL4.flank_3k.pi.combine", header = T)
tb <- read.table(args[1], header = T)

pdf(paste(args, "pdf", sep='.'), width=8, height = 5)

    #geom_line() +
#    geom_smooth(span = 0.001, se=F) +
#    geom_vline(xintercept = 4073880) +
#    geom_vline(xintercept = 4077946)  +
#    geom_smooth(span = 0.2, se=F) +
#    scale_y_continuous(expand = c(0,0), limits=c(0,0.4))+
    #geom_smooth(span = 0.1, se=F) +
min_val <- 4072002
max_val <- 4079753
ggplot(tb, aes(Loci, Pi,  color = Group)) +
    scale_x_continuous(breaks=seq(min_val, max_val, 800))+
    coord_cartesian(ylim=c(0,0.6)) +
    geom_smooth(span = 0.1, se=F) +
    scale_y_continuous(expand = c(0,0))+
    theme_bw()

dev.off()
