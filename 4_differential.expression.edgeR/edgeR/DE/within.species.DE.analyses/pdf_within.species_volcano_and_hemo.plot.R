library(tidyverse)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(LSD)
library(MASS)
library(viridis)
library(ggnewscale)
library(mdthemes)
library(VennDiagram)
library(ggtext)

rm(list = ls(all.names = T))



plot.colours<-c(
  "Cis only"="#000000",
  "Cis*"="#333333",
  "Cis**"="#707070",
  "Trans only"="#C7050F",
  "Trans*"="#FF5470",
  "Trans**"="#F88379",
  "Cis+Trans"="#8C12E2",
  "CisxTrans"="#2CA03D",
  "Compensatory"="#FFAB00",
  "Conserved"="gray90",
  "Ambiguous"="#0DD9E7",



  "Down regulated "="dodgerblue",
  "Down regulated"="dodgerblue",
  "Down regulated up LAM3 down PLASM1"="dodgerblue4",
  "Down regulated down LAM3 up PLASM1"="darkcyan",

  "Up regulated "="firebrick",
  "Up regulated"="firebrick",
  "Up regulated up LAM3 down PLASM1"="lightcoral",
  "Up regulated down LAM3 up PLASM1"="hotpink",

  "No significant difference "="black",
  "No significant difference"="black",
  "No significant difference up LAM3 down PLASM1"="gray45",
  "No significant difference down LAM3 up PLASM1"="gray"
)

load("Gene_info.RData")

DE.within <- bind_rows(read.csv("DE.result.fabo.response.within.SIM.SEC.HYB.csv"),
                       read.csv("DE.result.hemo.response.within.SIM.SEC.HYB.csv"))%>%
  mutate(species=factor(species,levels = c("*D. simulans*","Hybrid","*D. sechellia*")),
         tissue=factor(tissue,levels = c('Haemocytes','Fat Body')))%>%arrange(species,tissue)
DE.within$DE_trt.vs.ctrl[DE.within$FDR<0.05&DE.within$logFC>0] <- "Up regulated"
DE.within$DE_trt.vs.ctrl[DE.within$FDR<0.05&DE.within$logFC<0] <- "Down regulated"
DE.within<-merge(x = DE.within,y = Gene_info%>%dplyr::select(-dsim_gene),by = "dmel_gene",all.x = T)
DE.within$tissue<-gsub(pattern = "Haemo",replacement = "Hemo",x = DE.within$tissue)
DE.within <- DE.within%>%
  mutate(species=factor(species,levels = c("*D. simulans*","Hybrid","*D. sechellia*")),
         tissue=factor(tissue,levels = c('Hemocytes','Fat Body')))%>%arrange(species,tissue)

####################################
######### Within Species volcano ###
####################################
pdf(file = "Within.species.response.to.wasp.volcano.plot.pdf",width = 10,height = 10)
ggplot(data = DE.within,
mapping = aes(x = logFC,
              y = -log10(PValue),
              colour=DE_trt.vs.ctrl))+
  geom_point()+
  facet_grid(tissue~species)+
  scale_colour_manual(values = plot.colours)+
  labs(x="Log<sub>2</sub> Fold Change of Gene Expression Against Immune Challenge",
       y="-Log<sub>10</sub>(P value)")+
  xlim(-15,15)+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        legend.position = 'none',
        strip.background = element_rect(fill='NA',colour = 'black'),
        panel.grid = element_blank(),
        strip.text = element_markdown(size = 20),
        axis.title.x = element_markdown(size = 20),
        axis.title.y = element_markdown(size = 20),
        axis.text = element_text(size = 20))

dev.off()






################################################
########## violin plots compare express patterns of LAM3.marker+/- genes in control/wasp treatments, within each species
################################################

genes.hemo.immunne <- unique(DE.within$dmel_gene[DE.within$FDR<0.05&DE.within$tissue=='Hemocytes'])




hemo.plot <- ggplot(data = DE.within%>%filter(SC_sig_DE.cat!='',tissue=='Hemocytes',species!="Hybrid"),
                    mapping = aes(x = SC_sig_DE.cat,
                                  y = logFC,
                                  fill = SC_sig_DE.cat))+
  geom_abline(slope = 0,intercept = 0)+
  scale_fill_manual(values = c("up LAM3 down PLASM1"="lightcoral","down LAM3 up PLASM1"="darkcyan"),
                    labels =c("Up regulated","Down regulated"),
                    name = "Expression in Mature Lamellocytes",
                    guide = guide_legend(direction = "horizontal",title.position = "top"))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white")+
  facet_wrap(~species)+
  # ylim(-5,6)+
  labs(x = "",y="Log<sub>2</sub> Fold Change on Immune Challenge in Hemocytes")+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        strip.background = element_rect(fill='NA',colour='black'),
        panel.grid = element_blank(),
        strip.text = element_markdown(size = 12),
        legend.title = element_markdown(size=12),
        legend.position = 'top',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_markdown(size = 12),
        legend.text = element_markdown(size = 12))

save(hemo.plot,file = "hemo.plot.RData")

