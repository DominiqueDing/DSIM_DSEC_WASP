library(tidyverse)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(ggnewscale)
library(mdthemes)
library(ggtext)
library(ggsignif)
library(ggpubr)
library(lme4)

rm(list = ls(all.names = T))

load("Gene_info.RData")
genes.hemo.asde <- read.table("ASDE.analysed.genes.for.haemocytes.txt")
genes.hemo.immune <- read.table("immune.genes.hemo.txt")
#######################################
#### compare logfc between species ####
#######################################
DE.between.f1.parents <- read.csv("DE.result.hemo.response.betweenHSim.betweenHsec.csv")

data.DE.between.f1.parents <- bind_rows(
  DE.between.f1.parents%>%dplyr::select(dmel_gene,treatment,logFC.hyb_sim)%>%mutate(comp="hyb_sim")%>%dplyr::rename(logfc =logFC.hyb_sim),
  DE.between.f1.parents%>%dplyr::select(dmel_gene,treatment,logFC.hyb_sec)%>%mutate(comp="hyb_sec")%>%dplyr::rename(logfc =logFC.hyb_sec),
)
data.DE.between.f1.parents.1<-merge.data.frame(x = data.DE.between.f1.parents,y = Gene_info%>%dplyr::select(-dsim_gene)%>%filter(SC_sig_DE.cat!=""),by = "dmel_gene")
data.DE.between.f1.parents.1$comp<-gsub(pattern = "hyb_sim",replacement = "Hybrid/*D. simulans*",x = data.DE.between.f1.parents.1$comp)
data.DE.between.f1.parents.1$comp<-gsub(pattern = "hyb_sec",replacement = "Hybrid/*D. sechellia*",x = data.DE.between.f1.parents.1$comp)
pdf("Violin_DE.between.F1.parents_SC.and.immune.pdf",width = 6,height = 8)
ggplot(data = data.DE.between.f1.parents.1%>%filter(dmel_gene%in%genes.hemo.immune$V1,
                                                    treatment!="Response"),
       mapping = aes(x = SC_sig_DE.cat,
                     y = logfc,
                     fill = SC_sig_DE.cat))+
  geom_abline(slope = 0,intercept = 0)+
  scale_fill_manual(values = c("up LAM3 down PLASM1"="lightcoral","down LAM3 up PLASM1"="darkcyan"),
                    labels =c("Up regulated","Down regulated"),
                    guide = guide_legend(direction = "horizontal",title.position = "top"))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white")+
  facet_grid(comp~treatment)+
  # ylim(-5,6)+
  labs(x = "",y="Log<sub>2</sub> FC",
       fill="Expression in Mature Lamellocytes")+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        strip.background = element_rect(fill='NA',colour='black'),
        panel.grid = element_blank(),
        strip.text.x = element_markdown(size = 18),
        strip.text.y = element_markdown(size = 18),
        legend.title = element_markdown(size=18),
        legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_markdown(size = 20),
        legend.text = element_text(size = 18))
dev.off()
