library(tidyverse)
library(ggpubr)
library(ggtext)
library(ggsignif)
library(ggh4x)
library(extrafont)
font_import()
fonts()

rm(list = ls(all.names = T))
load("fabo.plot.RData")
load("hemo.plot.RData")

pdf("Figure7.pdf",width = 13,height = 8)
ggarrange(
  ggarrange(NULL,hemo.plot,NULL,nrow = 3,heights = c(0.05,1,0.05)),NULL,ggarrange(
    fabo.plot.1,
    NULL,
    ggarrange(
      plot.imd,plot.toll,plot.jak,
      ncol = 3, widths = c(3.6,3.6,2),
      labels = c("Imd","Toll","JAK-STAT"),hjust = c(-7,-9,-1.2),font.label = list(face='plain',size=12)
    ),
    nrow = 3,labels = c("B","C"),hjust = 1,heights = c(1,0.1,1),font.label = list(family="sans",face="plain")
  ),
  
  ncol = 3,labels = "A",hjust = 1,widths = c(1,0.1,3),font.label = list(family="sans",face="plain")
)+theme(plot.margin = margin(0.1,0.1,0.1,0.35, "cm")) 

dev.off()

