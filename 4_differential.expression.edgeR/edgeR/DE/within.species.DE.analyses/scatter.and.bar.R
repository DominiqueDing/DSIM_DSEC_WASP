library(tidyverse)
library(plyr)
library(dplyr)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(LSD)
library(MASS)
library(viridis)
library(ggnewscale)
library(mdthemes)
library(VennDiagram)
library(ggtext)
library(ggsignif)
library(ggpubr)
library(ggh4x)

rm(list = ls(all.names = T))



within.hemo <- read.csv("DE.result.hemo.response.within.SIM.SEC.HYB.csv")
within.fabo <- read.csv("DE.result.fabo.response.within.SIM.SEC.HYB.csv")
within.hemo$tissue<-gsub(pattern = "Haemo",replacement = "Hemo",x = within.hemo$tissue)
within.data <- bind_rows(within.hemo,within.fabo)%>%
  mutate(species=factor(species,levels = c("*D. simulans*","Hybrid","*D. sechellia*")),
         tissue=factor(tissue,levels = c('Hemocytes','Fat Body')))%>%arrange(species,tissue)

genes.hemo <- unique(within.hemo$dmel_gene[within.hemo$FDR<0.05])
genes.fabo <- unique(within.fabo$dmel_gene[within.fabo$FDR<0.05])
within.data.filter <- bind_rows(within.hemo%>%filter(dmel_gene%in%genes.hemo),within.fabo%>%filter(dmel_gene%in%genes.fabo))



log.plot<-ggplot(data = merge.data.frame(x = within.data.filter%>%filter(species=='*D. simulans*'),
                                         y = within.data.filter%>%filter(species=='*D. sechellia*'),
                                         by = c('dmel_gene','tissue'),
                                         suffixes = c('.sim','.sec'))%>%mutate(tissue=factor(tissue,levels = c('Hemocytes','Fat Body')),treatment='Response'),
                 mapping = aes(x = logFC.sim,y = logFC.sec))+
  geom_point(alpha=0.4)+
  geom_vline(xintercept = 0,linetype=3)+geom_hline(yintercept = 0,linetype=3)+geom_abline(slope = 1,intercept = 0,linetype=3)+
  xlim(-15,15)+ylim(-15,15)+
  facet_grid(rows=vars(tissue),cols = vars(treatment))+
  labs(x='Log<sub>2</sub> Fold Change of *D. simulans* <br>Expression on Immune Challenge',
       y='Log<sub>2</sub> Fold Change of *D. sechellia* Expression on Immune Challenge')+
  theme(strip.background = element_rect(fill = NA,colour = 'black'),
        strip.text = element_text(size=20),
        # strip.text.y = element_blank(),
        axis.text = element_text(size=20),
        panel.spacing = unit(x = 0,units = 'lines'),
        panel.border = element_rect(fill = 'NA',linetype = 1),
        panel.background = element_blank(),
        axis.title.x = element_markdown(size=20),axis.title.y = element_markdown(size=20))

pdf(file = 'DE.genes.responding.to.challenge.scatter.plot.pdf',width = 4.7,height = 10)
log.plot
dev.off()

cor.test(within.data.filter$logFC[within.data.filter$species=='*D. simulans*'&within.data.filter$tissue=='Hemocytes'],
         within.data.filter$logFC[within.data.filter$species=='*D. sechellia*'&within.data.filter$tissue=='Hemocytes'],method ="spearman")
cor.test(within.data.filter$logFC[within.data.filter$species=='*D. simulans*'&within.data.filter$tissue=='Fat Body'],
         within.data.filter$logFC[within.data.filter$species=='*D. sechellia*'&within.data.filter$tissue=='Fat Body'],method ="spearman")













#################################################
#################################################
#####   bar plot     ###########
#################################################
#################################################
reac.total.fabo <- read.csv("DE.result.fabo.response.betweenSIMSEC.csv")%>%filter(treatment=='Response')
reac.total.hemo <- read.csv("DE.result.hemo.response.betweenSIMSEC.csv")%>%filter(treatment=='Response')
within.data<-within.data%>%filter(species!='Hybrid')

### fabo ############

gene.plot.fabo <-merge.data.frame(x = within.data%>%filter(tissue=='Fat Body'),
                                  y = reac.total.fabo,
                                  by = c('dmel_gene','tissue'),suffixes = c("_trt.vs.ctrl",".total"),all.x = T)%>%
  mutate(species=factor(species,levels = c('*D. simulans*','Hybrid','*D. sechellia*')))%>%arrange(species)
for (i in c('within.data','gene.plot.fabo')) {
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "SPE.Dsim_GD21025",replacement = "SPE",x = .GlobalEnv[[i]]$dmel_gene)
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "AttA.AttB.Dsim_AttA",replacement = "AttA",x = .GlobalEnv[[i]]$dmel_gene)
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "AttA.AttB.Dsim_AttB",replacement = "AttB",x = .GlobalEnv[[i]]$dmel_gene)
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "CecA1.CecA2.CecC.Dsim_CecA1",replacement = "CecA1",x = .GlobalEnv[[i]]$dmel_gene)
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "CecA1.CecA2.CecC.Dsim_CecC",replacement = "CecC",x = .GlobalEnv[[i]]$dmel_gene)
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "TotM.TotF.Dsim_GD21729",replacement = "TotM",x = .GlobalEnv[[i]]$dmel_gene)
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "CG13749.CG33470.IMPPP.Dsim_GD10615",replacement = "BaraA2",x = .GlobalEnv[[i]]$dmel_gene)
  .GlobalEnv[[i]]$dmel_gene<-gsub(pattern = "TotA.TotB.Dsim_GD19384",replacement = "TotB",x = .GlobalEnv[[i]]$dmel_gene)
}
gene.plot.fabo$sig.within.sp[gene.plot.fabo$FDR_trt.vs.ctrl<0.05] <- '*'
gene.plot.fabo$sig.within.sp[gene.plot.fabo$FDR_trt.vs.ctrl<0.01] <- '**'
gene.plot.fabo$sig.within.sp[gene.plot.fabo$FDR_trt.vs.ctrl<0.001] <- '***'
gene.plot.fabo$sig.total[gene.plot.fabo$FDR.total<0.05] <- '*'
gene.plot.fabo$sig.total[gene.plot.fabo$FDR.total<0.01] <- '**'
gene.plot.fabo$sig.total[gene.plot.fabo$FDR.total<0.001] <- '***'


gene.fabo.immune.func <- c("Dro","Mtk","Def","AttA","CecC","BaraA2",
                           "BomBc1","BomS1","BomS4",
                           "TotM","TotF","TotA","TotB","TotX",
                           # "SPE","CG30098",
                           "Sp212",
                           #"PPO3",
                           "edin")
data = gene.plot.fabo%>%filter(dmel_gene%in%gene.fabo.immune.func,species!='Hybrid')%>%mutate(dmel_gene=factor(dmel_gene,levels = c("Dro","Mtk","Def","AttA","CecC","BaraA2",
                                                                                                                                    "BomBc1","BomS1","BomS4",
                                                                                                                                    "TotM","TotF","TotA","TotB","TotX",
                                                                                                                                    # "SPE","CG30098",
                                                                                                                                    "Sp212",
                                                                                                                                    #"PPO3",
                                                                                                                                    "edin")))%>%arrange(dmel_gene)

fabo.plot.1 <- ggplot(data = data,
                      mapping = aes(x = dmel_gene,y = logFC_trt.vs.ctrl,fill=species))+
  geom_col(position = position_dodge(),width = 0.7)+
  geom_text(aes(label=sig.within.sp),position = position_dodge(width=0.7),vjust = -0.1)+
  geom_signif(annotations = c("*"),
              y_position = c(11.5), xmin=c(9.8), xmax=c(10.2))+
  annotate(geom = 'segment',x = 1,xend = 6,y = 15,yend = 15)+
  annotate(geom = 'text',x = 3.5,y = 16,label='AMP genes',size=5)+
  annotate(geom = 'segment',x = 7,xend = 9,y = 15,yend = 15)+
  annotate(geom = 'text',x = 8,y = 16,label='Bomanin genes',size=5)+
  annotate(geom = 'segment',x = 10,xend = 12,y = 15,yend = 15)+
  annotate(geom = 'text',x = 11,y = 16,label='Turandot genes',size=5)+
  annotate(geom = 'segment',x = 12.7,xend = 13.3,y = 15,yend = 15)+
  annotate(geom = 'text',x = 13,y = 16,label='edin',size=5)+
  expand_limits(y=c(0,13))+
  labs(y='Log<sub>2</sub> Fold Change on Immune Challenge')+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.y = element_markdown(size=12),
        axis.text.y = element_text(size=12),
        legend.position='top',legend.justification='right',legend.direction='horizontal',
        legend.title = element_blank(),
        legend.text = element_markdown(size=12))




core.imd <- c('imd','Fadd','Dredd','key','Rel')
core.toll <- c('SPE','spz','Tl','Dif')
core.jak <- c('hop','Stat92E')
data <- bind_rows(
  gene.plot.fabo%>%filter(dmel_gene%in%core.imd)%>%mutate(dmel_gene=factor(dmel_gene,levels =  core.imd),path='Imd Signaling Pathway'),
  gene.plot.fabo%>%filter(dmel_gene%in%core.toll)%>%mutate(dmel_gene=factor(dmel_gene,levels =  core.toll),path='Toll Signaling Pathway'),
  gene.plot.fabo%>%filter(dmel_gene%in%core.jak)%>%mutate(dmel_gene=factor(dmel_gene,levels =  core.jak),path='JAK-STAT Signaling Pathway')
)%>%filter(species!='Hybrid')

plot.imd<-ggplot(data = data%>%filter(path=='Imd Signaling Pathway'),
                 mapping = aes(x = dmel_gene,y = logFC_trt.vs.ctrl,fill=species))+
  geom_col(position = position_dodge(),width = 0.7)+
  geom_text(aes(label=sig.within.sp),position = position_dodge(width=0.7),vjust = -0.1)+
  geom_signif(annotations = c("*","*","*"),
              y_position = c(0.8,1.3,1.7), xmin=c(0.8,2.8,3.8), xmax=c(1.2,3.2,4.2))+
  expand_limits(y=c(0,2))+
  labs(y='Log<sub>2</sub> Fold Change on Immune Challenge')+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(),
        axis.title.y = element_markdown(size=12),
        axis.text = element_text(size=12),
        legend.position = 'none',
        # legend.position='top',legend.justification='right',legend.direction='horizontal',
        legend.title = element_blank(),
        legend.text = element_markdown())
plot.toll<-ggplot(data = data%>%filter(path=='Toll Signaling Pathway'),
                 mapping = aes(x = dmel_gene,y = logFC_trt.vs.ctrl,fill=species))+
  geom_col(position = position_dodge(),width = 0.7)+
  geom_text(aes(label=sig.within.sp),position = position_dodge(width=0.7),vjust = -0.1)+
  geom_signif(annotations = c("*"),
              y_position = c(2.3), xmin=c(0.8), xmax=c(1.2))+
  expand_limits(y=c(0,2.5))+
  # labs(y='Log<sub>2</sub> Fold Change on Immune Challenge')+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        legend.position = 'none',
        # legend.position='top',legend.justification='right',legend.direction='horizontal',
        legend.title = element_blank(),
        legend.text = element_markdown(size=12))
plot.jak<-ggplot(data = data%>%filter(path=='JAK-STAT Signaling Pathway'),
                  mapping = aes(x = dmel_gene,y = logFC_trt.vs.ctrl,fill=species))+
  geom_col(position = position_dodge(),width = 0.7)+
  geom_text(aes(label=sig.within.sp),position = position_dodge(width=0.7),vjust = -0.1)+
  # geom_signif(annotations = c("*", "*","*","*","**","**"),
  #             y_position = c(1.2,1.65,0.7,2,2.2,1.6), xmin=c(0.8,1.8,2.8,3.8,4.8,5.8), xmax=c(1.2,2.2,3.2,4.2,5.2,6.2))+
  expand_limits(y=c(0,1.5))+
  # labs(y='Log<sub>2</sub> Fold Change on Immune Challenge')+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        legend.position = 'none',
        # legend.position='top',legend.justification='right',legend.direction='horizontal',
        legend.title = element_blank(),
        legend.text = element_markdown(size=12))



save(fabo.plot.1, plot.imd,plot.jak,plot.toll,file = "fabo.plot.RData")



  
  