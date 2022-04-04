
library(tidyverse)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(ggtext)
library(markdown)
library(viridis)
library(hrbrthemes)
library(scales)



remove(list = ls(all.names = T))


plot.colours<-c(
  "*Cis*-diverged"="#000000",
  "*Trans*-diverged"="#C7050F",
  # "*Cis*+*Trans*"="#8C12E2",
  # "*CisxTrans*"="#2CA03D",
  # "Compensatory"="#FFAB00",
  "Both *cis*- and *trans*-diverged"="#F2D77D",
  "Ambiguous"="#0DD9E7",
  "Conserved"="#ACBAC3"
)








ana.ctrl <- read.csv("ASDE.result.hemo.ctrl.csv")
ana.wasp <- read.csv("ASDE.result.hemo.wasp.csv")
ana.reac <- read.csv("ASDE.result.hemo.reac.csv")
DE.within <- read.csv("DE.result.hemo.response.within.SIM.SEC.HYB.csv")

genes.immune.responsive <- unique(DE.within$dmel_gene[DE.within$FDR<0.05])


for (i in c('ana.ctrl','ana.wasp','ana.reac')) {
  .GlobalEnv[[i]] <- .GlobalEnv[[i]]%>%filter(.GlobalEnv[[i]]$dmel_gene%in%genes.immune.responsive)
  
  .GlobalEnv[[paste0(i,".2.1")]] <- .GlobalEnv[[i]]
  
  .GlobalEnv[[paste0(i,".2.1")]]$fill.2 <- .GlobalEnv[[paste0(i,".2.1")]]$Cat.main
  .GlobalEnv[[paste0(i,".2.1")]]$fill.2 <- gsub(pattern = "Cis related",replacement = "*Cis*-diverged",x = .GlobalEnv[[paste0(i,".2.1")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.1")]]$fill.2 <- gsub(pattern = "Trans related",replacement = "*Trans*-diverged",x = .GlobalEnv[[paste0(i,".2.1")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.1")]]$fill.2 <- gsub(pattern = "Cis\\+Trans",replacement = "Both *cis*- and *trans*-diverged",x = .GlobalEnv[[paste0(i,".2.1")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.1")]]$fill.2 <- gsub(pattern = "CisxTrans",replacement = "Both *cis*- and *trans*-diverged",x = .GlobalEnv[[paste0(i,".2.1")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.1")]]$fill.2 <- gsub(pattern = "Compensatory",replacement = "Both *cis*- and *trans*-diverged",x = .GlobalEnv[[paste0(i,".2.1")]]$fill.2)
  
  .GlobalEnv[[paste0(i,".2.1")]]$Cat.main <- gsub(pattern = "Cis related",replacement = "*Cis*",x = .GlobalEnv[[paste0(i,".2.1")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.1")]]$Cat.main <- gsub(pattern = "Trans related",replacement = "*Trans*",x = .GlobalEnv[[paste0(i,".2.1")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.1")]]$Cat.main <- gsub(pattern = "Cis\\+Trans",replacement = "*Cis*",x = .GlobalEnv[[paste0(i,".2.1")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.1")]]$Cat.main <- gsub(pattern = "CisxTrans",replacement = "*Cis*",x = .GlobalEnv[[paste0(i,".2.1")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.1")]]$Cat.main <- gsub(pattern = "Compensatory",replacement = "*Cis*",x = .GlobalEnv[[paste0(i,".2.1")]]$Cat.main)
  
  .GlobalEnv[[paste0(i,".2.2")]] <- .GlobalEnv[[i]]
  
  .GlobalEnv[[paste0(i,".2.2")]]$fill.2 <- .GlobalEnv[[paste0(i,".2.2")]]$Cat.main
  .GlobalEnv[[paste0(i,".2.2")]]$fill.2 <- gsub(pattern = "Cis related",replacement = "*Cis*-diverged",x = .GlobalEnv[[paste0(i,".2.2")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.2")]]$fill.2 <- gsub(pattern = "Trans related",replacement = "*Trans*-diverged",x = .GlobalEnv[[paste0(i,".2.2")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.2")]]$fill.2 <- gsub(pattern = "Cis\\+Trans",replacement = "Both *cis*- and *trans*-diverged",x = .GlobalEnv[[paste0(i,".2.2")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.2")]]$fill.2 <- gsub(pattern = "CisxTrans",replacement = "Both *cis*- and *trans*-diverged",x = .GlobalEnv[[paste0(i,".2.2")]]$fill.2)
  .GlobalEnv[[paste0(i,".2.2")]]$fill.2 <- gsub(pattern = "Compensatory",replacement = "Both *cis*- and *trans*-diverged",x = .GlobalEnv[[paste0(i,".2.2")]]$fill.2)
  
  .GlobalEnv[[paste0(i,".2.2")]]$Cat.main <- gsub(pattern = "Cis related",replacement = "*Cis*",x = .GlobalEnv[[paste0(i,".2.2")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.2")]]$Cat.main <- gsub(pattern = "Trans related",replacement = "*Trans*",x = .GlobalEnv[[paste0(i,".2.2")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.2")]]$Cat.main <- gsub(pattern = "Cis\\+Trans",replacement = "*Trans*",x = .GlobalEnv[[paste0(i,".2.2")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.2")]]$Cat.main <- gsub(pattern = "CisxTrans",replacement = "*Trans*",x = .GlobalEnv[[paste0(i,".2.2")]]$Cat.main)
  .GlobalEnv[[paste0(i,".2.2")]]$Cat.main <- gsub(pattern = "Compensatory",replacement = "*Trans*",x = .GlobalEnv[[paste0(i,".2.2")]]$Cat.main)
  
  .GlobalEnv[[paste0(i,".2")]] <- bind_rows(.GlobalEnv[[paste0(i,".2.1")]],.GlobalEnv[[paste0(i,".2.2")]])%>%distinct()
}
ana.con.2 <- merge.data.frame(x = ana.wasp.2%>%select(-fill.2),y = ana.ctrl.2%>%select(dmel_gene,fill.2)%>%distinct(),by = "dmel_gene",all.x  = T)

data = bind_rows(ana.ctrl.2%>%mutate(treatment="Control"),
                 ana.wasp.2%>%mutate(treatment="Wasp"),
                 ana.reac.2%>%mutate(treatment="Response (wasp - control)"),
                 ana.con.2%>%mutate(treatment="con"),
)%>%
  mutate(treatment = factor(treatment,levels = c("Control","Wasp","Response (wasp - control)","con")))%>%
  mutate(Cat.main=factor(Cat.main,
                         levels = c("Conserved",
                                    "Ambiguous",
                                    "Compensatory",
                                    "*CisxTrans*",
                                    "*Cis*+*Trans*",
                                    "*Trans*",
                                    "*Cis*")),
         fill.2=factor(fill.2,levels = c("Conserved","Ambiguous","Both *cis*- and *trans*-diverged","*Cis*-diverged","*Trans*-diverged")))%>%arrange(fill.2,Cat.main)


pdf("cis.trans.log.plot.hemo_immune.responsive.genes.pdf",width = 15,height = 6.3)
ggplot(data = data%>%filter(treatment!='con'),
       mapping = aes(x = logFC.total,y = logFC.cis))+
  scale_colour_manual(values = plot.colours)+
  geom_abline(slope = 1,intercept = 0,linetype = 3)+geom_abline(slope = 0,intercept = 0,linetype = 3)+geom_vline(xintercept = 0,linetype = 3)+
  annotate(geom = "text",x = 14.5,y = 0.5,label = "italic(trans)~only",parse=T,colour = "red",fontface = "bold", size = 5)+
  annotate(geom = "text",x = 14.5,y = 16,label = "italic(cis)~only",parse=T,colour = "black", fontface = "bold",size = 5, angle = 45)+
  geom_point(aes(colour = fill.2))+
  labs(
    x="Log<sub>2</sub>(Parent<sub>*sim*</sub>/Parent<sub>*sec*</sub>)",
    y="Log<sub>2</sub>(Hybrid<sub>*sim*</sub>/Hybrid<sub>*sec*</sub>)")+
  theme_bw()+
  facet_grid(cols = vars(treatment))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size = 20),
        legend.position = "none",
        panel.grid = element_blank())+
  xlim(-17,17)+ylim(-17,17)
dev.off()

pdf("cis.trans.bar.plot.hemo_immune.responsive.genes.pdf",width = 15,height = 2.5)
ggplot(data = data%>%filter(treatment!='con',Cat.main!='Conserved'),
       mapping = aes())+
  scale_fill_manual(values = plot.colours)+
  geom_bar(aes(y = Cat.main, fill = fill.2),position = position_stack(reverse = T))+
  facet_grid(cols = vars(treatment))+
  labs(
    x="Number of genes",
    y="")+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text.y = element_markdown(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(1.2,0.5,0.5,0.5),"cm"),
        legend.position=c(0.55, 1.2),
        legend.title = element_blank(),
        # legend.position = "none",
        panel.grid = element_blank())
dev.off()




data.hist<-bind_rows(ana.ctrl.2.1%>%mutate(treatment='Control'),
                     ana.wasp.2.1%>%mutate(treatment='Wasp'),
                     ana.reac.2.1%>%mutate(treatment='Response'))%>%
  filter(Cat.main=="*Cis*"|Cat.main=="*Trans*")%>%
  mutate(treatment=factor(treatment,levels = c('Control','Wasp','Response')))
data.space <- data.frame(dmel_gene=rep(x = 'space',90),
                        logFC.total=rep(seq(0,14),6),
                        Cat.main=rep(c(rep('*Cis*',15),rep('*Trans*',15)),3),
                        treatment=c(rep('Control',30),rep('Wasp',30),rep('Response',30)))
data.hist<-merge.data.frame(x = data.hist,y = data.space,by = names(data.space),all = T)

pdf('hemo.cis.trans.number.genes.hist.plot_immune.responsive.genes.pdf',width = 15,height = 5)
ggplot(data = data.hist%>%mutate(logFC.total=abs(logFC.total)))+
  geom_histogram(aes(x = logFC.total,fill = Cat.main),colour='black',position = position_dodge(),bins = 16)+
  scale_fill_manual(values = list("*Cis*" = "#000000","*Trans*" = "#C7050F"))+
  facet_grid(cols = vars(treatment))+
  scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
                labels = trans_format(trans = "log10",format = math_format(.x)),
                limits=c(NA,1000),
                oob=squish_infinite)+
  labs(y='Log<sub>10</sub> Number of Genes',
       x='Magnitude of expression divergence (log<sub>2</sub>Parent<sub>*sim*</sub>/Parent<sub>*sec*</sub>)')+
  theme_bw()+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0,'lines'),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position=c(0.92, 0.9),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size=20),
        strip.background = element_rect(fill=F,colour='black'),
        strip.text.x = element_markdown(size=20))
dev.off()
  
  









genes.ctrl.cis <- unique(ana.ctrl.2$dmel_gene[ana.ctrl.2$Cat.main=='*Cis*'])
genes.ctrl.only.cis <- unique(ana.ctrl.2$dmel_gene[ana.ctrl.2$fill.2=='*Cis*-diverged'])
genes.wasp.cis <- unique(ana.wasp.2$dmel_gene[ana.wasp.2$Cat.main=='*Cis*'])
genes.wasp.only.cis <- unique(ana.wasp.2$dmel_gene[ana.wasp.2$fill.2=='*Cis*-diverged'])

genes.ctrl.trans <- unique(ana.ctrl.2$dmel_gene[ana.ctrl.2$Cat.main=='*Trans*'])
genes.ctrl.only.trans <- unique(ana.ctrl.2$dmel_gene[ana.ctrl.2$fill.2=='*Trans*-diverged'])
genes.wasp.trans <- unique(ana.wasp.2$dmel_gene[ana.wasp.2$Cat.main=='*Trans*'])
genes.wasp.only.trans <- unique(ana.wasp.2$dmel_gene[ana.wasp.2$fill.2=='*Trans*-diverged'])


cor.test(ana.ctrl$logFC.cis,ana.ctrl$logFC.total)
cor.test(ana.ctrl$logFC.trans,ana.ctrl$logFC.total)

cor.test(ana.wasp$logFC.cis,ana.wasp$logFC.total)
cor.test(ana.wasp$logFC.trans,ana.wasp$logFC.total)

cor.test(ana.reac$logFC.cis,ana.reac$logFC.total)
cor.test(ana.reac$logFC.trans,ana.reac$logFC.total)



