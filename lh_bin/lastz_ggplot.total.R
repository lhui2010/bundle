#!/usr/bin/env Rscript

library("ggplot2")
#install.packages("scales")
library("scales")

#---
theme_custom <- function() {
  theme_bw() + # note ggplot2 theme is used as a basis
  theme(plot.title = element_text(size = 14, face = "bold",
                                  hjust = .5,
                                  margin = margin(t = 5, b = 15)),
        plot.caption = element_text(size = 11, hjust = 0,
                                    margin = margin(t = 15)),
        panel.grid.major = element_line(colour = "grey88"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10),
                                    size = 13, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10),
                                    size = 13, face = "bold"))
}
theme_Publication <- function(base_size=14, base_family="sans") {
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
               #jlegend.position = "right",
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
#---

args = commandArgs(trailingOnly=TRUE)
output = paste(args[1], ".gg.pdf", sep="")
#pdf(paste(args[1], ".gg.pdf", sep=""))

fname = args[1]

test<-read.table(fname, header=TRUE,comment.char = "")

forward <- subset(test, test$strand2 == '+')
backward <- subset(test, test$strand2 == '-')

neworder <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")
library(plyr)  ## or dplyr (transform -> mutate)
test2 <- arrange(transform(test,
             chr=factor(chr,levels=neworder)),chr)

b <- ggplot(test2, aes(x = zstart1, y = zstart2.))

b + geom_segment(aes(xend=end1, yend=end2., color = strand2), size=1, show.legend = FALSE) +
 scale_x_continuous(breaks = breaks_width(200000000), labels = label_number(scale = 1/1000000, suffix="Mb")) +
 scale_y_continuous(breaks = breaks_width(50000000), labels = label_number(scale = 1/1000000, suffix="Mb")) +
 labs(x = "A188", y = "B73") + 
 #labs(x = "A188", y = "B73") + 
 #facet_wrap(~ chr) +
 facet_grid(cols = vars(chr)) +
 theme_Publication() 
 #theme_custom()

aspect_ratio <- 10
height <- 7
ggsave(output , height = 2 , width = 2 * aspect_ratio)


##plot(NA, xlim=c(0,max(test$end1)), ylim=c(0,max(test$end2.)),
##     xlab="ref", ylab="qry", main=fname)
##
##segments(forward$zstart1,forward$zstart2.,forward$end1, forward$end2., col='blue', lwd=1) 
##segments(backward$zstart1,backward$zstart2.,backward$end1, backward$end2., col='red', lwd=1 )
#dev.off()
