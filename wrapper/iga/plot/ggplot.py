"""
ggplot wrappers
"""
from iga.apps.base import emain, sh, rscript

theme_publication_r = r"""
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
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506", 
      "#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506",
      "#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
"""
# 0 input_table
# 1 x=xx
# 2 y=xx
# 3 fill=xx
# 4 etc string like theme and others
barplot_r = r"""
library(ggplot2);
a=read.table("{0}", header=T, row.names=NULL); 
ggplot(a, aes(x={1}, y={2}, fill={3})) + geom_bar(stat="identity", position=position_dodge()) {4} 
ggsave("{0}.pdf", width = 5, height = 5)
"""


def barplot(table=None, x='', y='', group='', theme='Publication', horizonal='F'):
    r"""
    barplot with ggplot2
    :param table: Input table like
                Group   Number  Regulation
                Root-Nodule     1654    Down
                Root-Nodule     1496    Up
                Root-Leaf       4099    Down
    :param x: default is the 1st colomn, can be specified by header
    :param y: default is the 2nd colomn, can be specified by header
    :param group: default is the 3rd colomn, can be specified by header
    :param theme: available themes(Publication, minimal)
    :param horizonal: whether to plot horizonally. (T|F default F)
    :return:
    """
    if x == '':
        with open(table) as fh:
            header = fh.readline()
            (x, y, group) = header.rstrip().split()
    etc = ''
    if theme != "":
        etc = "+theme_" + theme + "()"
    if horizonal == 'T':
        etc += '+ coord_flip()'
    cmd = theme_publication_r + barplot_r.format(table, x, y, group, etc)
    rscript(cmd)


if __name__ == "__main__":
    emain()
