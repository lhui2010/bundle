#!/usr/bin/env Rscript
###Input example
#Group	Loci	Fst
#Aus_vs_Indica	447	0.028819868067541093
#Aus_vs_Indica	100933	0.038394295853923385
#Aus_vs_Indica	200973	0.042556573585786975
#Aus_vs_Indica	301162	0.03157124525044273
#Aus_vs_Indica	401213	0.022624335507077185
#Aus_vs_Indica	501215	0.03651418067888159
#Aus_vs_Indica	601288	0.03289836343059585
#Aus_vs_Indica	701306	0.04272873051736839
#Aus_vs_Indica	801357	0.04510758247521361



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

#    geom_smooth(span = 0.2, se=F) +
#    coord_cartesian(ylim=c(0,0.4)) +
#    scale_y_continuous(expand = c(0,0))+
#    scale_x_continuous(breaks=seq(1, 35935812, 1000000))+
#    scale_y_continuous(expand = c(0,0))+

#library("sitools")
#    scale_x_log10("Message Size [Byte]", labels=f2si) +
ggplot(tb, aes(Loci, Fst)) +
    geom_line() +
    geom_smooth(span = 0.001, se=F) +
    geom_segment(aes(x=4073880, xend=4073880, y=0.25, yend=0.23), size = 0.5,
            arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text", x = 4073880, y = 0.27, label = "OsSPL4") +
    facet_wrap( ~ Group) +
    theme_bw()

dev.off()
