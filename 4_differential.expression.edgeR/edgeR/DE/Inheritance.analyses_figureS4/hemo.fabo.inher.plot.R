library(tidyverse)
library(ggrepel)
library(ggtext)
library(markdown)
library(ggpubr)

rm(list = ls(all.names = T))

plot.colours<-c(
  '*D. simulans* dominant' = "#F8766D",
  '*D. sechellia* dominant' = "#00BFC4",
  "Additive" = "#FFC300",
  "Overdominant" = "#089911",
  "Underdominant" = "#9B5DE5",
  "Conserved" = "#778DA9"
  
)
plot.colours.2<-c(
  '*D. simulans*' = "#F8766D",
  '*D. sechellia*' = "#00BFC4",
  'Hybrid' = '#089911'
  
)


data.hemo <- read.csv("inheritance.DE.between.species_hemo.csv")
data.fabo <- read.csv("inheritance.DE.between.species_fabo.csv")
genes.immune.hemo <- read.table("immune.genes.hemo.txt")
genes.immune.fabo <- read.table("immune.genes.fabo.txt")

data <- bind_rows(data.hemo%>%filter(dmel_gene%in%genes.immune.hemo$V1), data.fabo%>%filter(dmel_gene%in%genes.immune.fabo$V1))
data$tissue<-gsub(pattern = "Haemocytes",replacement = "Hemocytes",x = data$tissue)






pdf(file = "inheritance.bar.plot.pdf",width = 7,height = 1.5)
ggplot(data = data%>%filter(inher!="Conserved",treatment=="Response")%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. simulans* dominant",
                                              "*D. sechellia* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wasp','Response')),
                tissue=factor(tissue,levels = c("Hemocytes","Fat Body")))%>%arrange(inher,treatment))+
  scale_fill_manual(values = plot.colours)+
  geom_bar(mapping = aes(y=inher,fill = inher))+
  facet_wrap(~tissue,scales = "free_x")+
  labs(x = 'Number of genes',y='')+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.text.y = element_markdown())
dev.off()





bar.hemo <- ggplot(data = data%>%filter(inher!="Conserved",tissue=="Hemocytes")%>%
                     mutate(inher=factor(inher,levels = c("Conserved",
                                                          "*D. simulans* dominant",
                                                          "*D. sechellia* dominant",
                                                          "Additive",
                                                          "Overdominant",
                                                          "Underdominant")),
                            treatment=factor(treatment,levels = c('Control','Wasp','Response')),
                            tissue=factor(tissue,levels = c("Hemocytes","Fat Body")))%>%arrange(inher,treatment))+
  scale_fill_manual(values = plot.colours)+
  geom_bar(mapping = aes(y=inher,fill = inher))+
  facet_grid(cols = vars(treatment),rows = vars(tissue))+
  labs(y = '')+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.text.y = element_markdown(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1.2))
bar.fabo <- ggplot(data = data%>%filter(inher!="Conserved",tissue=="Fat Body")%>%
                     mutate(inher=factor(inher,levels = c("Conserved",
                                                          "*D. simulans* dominant",
                                                          "*D. sechellia* dominant",
                                                          "Additive",
                                                          "Overdominant",
                                                          "Underdominant")),
                            treatment=factor(treatment,levels = c('Control','Wasp','Response')),
                            tissue=factor(tissue,levels = c("Hemocytes","Fat Body")))%>%arrange(inher,treatment))+
  scale_fill_manual(values = plot.colours)+
  geom_bar(mapping = aes(y=inher,fill = inher))+
  facet_grid(cols = vars(treatment),rows = vars(tissue))+
  labs(y = '', x = 'Number of genes')+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        strip.text.x = element_blank(),
        axis.text.y = element_markdown(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1.2))

pdf(file = "supp_inheritance.bar.plot.pdf",width = 5,height = 3)
ggarrange(bar.hemo,bar.fabo,ncol = 1,nrow = 2)
dev.off()











pdf(file = "inheritance.bar.plot.pdf",width = 7,height = 3)

dev.off()



pdf(file = "inheritance.scatter.plot.pdf",width = 5,height = 3.5)
ggplot(data = data%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. simulans* dominant",
                                              "*D. sechellia* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wasp','Response')),
                tissue=factor(tissue,levels = c("Hemocytes","Fat Body")))%>%arrange(inher,treatment,tissue),                  
       mapping = aes(x = `logFC.hyb_sim`,y = `logFC.hyb_sec`,colour = inher))+
  scale_colour_manual(values = plot.colours)+
  geom_point(size = 1,shape=1)+
  geom_vline(xintercept = 0,colour='light grey',linetype=3)+geom_hline(yintercept = 0,colour='light grey',linetype=3)+
  xlim(-12,19)+ylim(-12,19)+
  guides(colour = guide_legend(override.aes =  list(size = 3.5,shape = 15)))+
  facet_grid(cols = vars(treatment),rows = vars(tissue),scales = "free_x")+
  labs(x = "Log<sub>2</sub>( Hybrid / *D. simulans* )",
       y = "Log<sub>2</sub>( Hybrid / *D. sechellia* )")+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())
dev.off()






