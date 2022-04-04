library(tidyverse)
library(ggtext)
library(ggrepel)
require(forcats)




remove(list = ls(all.names = T))
# load the RNA sequencing read counts for constitutive and conserved exon gene counts:
load("RNAseq.gene.counts.raw.RData")

colour <- list(
  '*D. simulans*' = "#f8766d",
  '*D. sechellia*' = "#00BFC4",
  'Hybrid' = "#089911"
  )

#######################################################
### filter out genes with>0.1 mis-mapped rates ########
#######################################################
## create genewise counts and filter out genes with>0.1 mis-mapped rates:
counts.genewise <- counts%>%
  select(species,counts,allele.species.consistency,dmel_gene)%>%
  group_by(species,allele.species.consistency,dmel_gene)%>%summarise(counts=sum(counts))
counts.genewise <- merge(x = counts.genewise[counts.genewise$allele.species.consistency=="Consistent with genotype",]%>%
                           select(species,dmel_gene,counts)%>%
                           mutate(counts.true=counts)%>%
                           select(-counts),
                         y = counts.genewise[counts.genewise$allele.species.consistency=="Inconsistent with genotype",]%>%
                           select(species,dmel_gene,counts)%>%
                           mutate(counts.false=counts)%>%
                           select(-counts),
                         by = c("species","dmel_gene"),all = T)%>%
  select(dmel_gene,species,counts.true,counts.false)
counts.genewise$f.p <- round(counts.genewise$counts.false/(counts.genewise$counts.false+counts.genewise$counts.true),3)

## get genes with high mis-mapped rates:
mismap_0.1_filtered_out_Genes<-unique(counts.genewise$dmel_gene[counts.genewise$f.p>0.1])

## filter:
counts <- counts[!counts$dmel_gene%in%mismap_0.1_filtered_out_Genes,]

## create libwise counts:
counts.libwise <- counts%>%
  select("lib","species","counts","allele.species.consistency")%>%
  group_by(lib,species,allele.species.consistency)%>%
  summarise(counts=sum(counts))
counts.libwise <- merge(
  x = counts.libwise[counts.libwise$allele.species.consistency=="Consistent with genotype",]%>%
    select(lib,species,counts)%>%mutate(counts.true = counts)%>%select(-counts),
  y = counts.libwise[counts.libwise$allele.species.consistency=="Inconsistent with genotype",]%>%
    select(lib,species,counts)%>%mutate(counts.false = counts)%>%select(-counts),
  by = c("lib","species"),all = T
)
counts.libwise[is.na(counts.libwise)]<-0
counts.libwise$f.p <- round(x = counts.libwise$counts.false/(counts.libwise$counts.false+counts.libwise$counts.true),
                            digits = 3)
counts.libwise <- bind_rows(
  counts.libwise%>%
    select(lib,species,counts.true,f.p)%>%
    mutate(counts=counts.true,allele.species.consistency="Consistent with genotype")%>%
    select(-counts.true),
  counts.libwise%>%
    select(lib,species,counts.false,f.p)%>%
    mutate(counts=counts.false,allele.species.consistency="Inconsistent with genotype")%>%
    select(-counts.false)
)

###adding lib info:
lib.info<-read.csv("library_information.csv")
counts.libwise.2 <- merge(x = counts.libwise,y = lib.info%>%select(lib,tissue),by = 'lib')%>%
  mutate(lib.label = paste0("<span style = 'color: ",ifelse(test = tissue=='Fat Body',yes = "#E0B000",no = "#D5202F"),";'>",lib,"</span>"))
counts.libwise.2.x<-counts.libwise.2%>%group_by(lib)%>%summarise(t.counts=sum(counts))
counts.libwise.2<-merge(x=counts.libwise.2,y=counts.libwise.2.x,all = T)%>%mutate(order=f.p*10e10+t.counts)
 

pdf("Allelic.check.1.2_Fig.S3.A.pdf",width = 13,height = 5)
ggplot(data = counts.libwise.2[counts.libwise.2$species!="Hybrid",])+
  geom_col(
    aes(x = reorder(lib.label,-order),
        y = counts,
        fill = allele.species.consistency),
    position = position_stack(reverse = T)
  )+
  scale_fill_manual(values = list("Inconsistent with genotype" = "grey60",
                                  "Consistent with genotype" = "grey25"))+
  geom_text(mapping = aes(x = lib.label,
                          y = counts+500000,
                          label=ifelse(test = (f.p>0&counts>4000000),yes = paste0(f.p*100,'%'),no = '')),
            size=4.5)+
  coord_flip()+
  facet_wrap(facets = ~species,scales = 'free_y')+
  theme_bw()+
  labs(title = "<span style = 'color: #E0B000;'>Fat Body library</span><br><span style = 'color: #D5202F;'>Haemocytes library</span>")+
  xlab("Library")+
  ylab("RNAseq read counts")+
  theme(strip.text = element_markdown(size = 24),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_markdown(size = 14),
        axis.text.y = element_markdown(size = 14),
        plot.title = element_markdown(size=20,hjust = 1),
        legend.position = 'none')+
  scale_y_continuous(breaks = c(0,2.5e+06,5.0e+06,7.5e+06,1.0e+07),
                     labels = c("0","2.5e+06","5.0e+06","7.5e+06","1.0e+07"))
dev.off()

counts.libwise.3 <- counts%>%
  mutate(allele=substr(x = allele,start = 1,stop = 4))%>%
  select(lib,allele,species,counts)%>%
  group_by(lib,allele,species)%>%summarise(counts=sum(counts))
counts.libwise.3 <- merge(x = counts.libwise.3[counts.libwise.3$allele=="dsim",]%>%
                       select(lib,species,counts)%>%mutate(counts.dsim = counts),
                     y = counts.libwise.3[counts.libwise.3$allele=="dsec",]%>%
                       select(lib,species,counts)%>%mutate(counts.dsec=counts),
                     by = c("lib","species"))%>%
  select(lib,species,counts.dsim,counts.dsec)%>%
  mutate(species = factor(x = species,levels = c("*D. simulans*","Hybrid","*D. sechellia*")))
pdf("Allelic.check.3.pdf_Fig.S3.B.pdf",width = 3.3,height = 3.3)
ggplot(data = counts.libwise.3)+
  scale_colour_manual(values = colour)+
  geom_point(mapping = aes(x = log2(counts.dsim),
                           y = log2(counts.dsec),
                           colour = species),size=1.5)+
  geom_abline(slope = 1,intercept = 0,linetype=3)+
  geom_text_repel(mapping = aes(x = log2(counts.dsim),
                                y = log2(counts.dsec),
                                label = ifelse(test = lib%in%c('B7','D7'),yes = as.character(lib),no = '')))+
  xlab("Log<sub>2</sub>(*D. simulans* allelic read counts)")+
  ylab("Log<sub>2</sub>(*D. sechellia* allelic read counts)")+
  # labs(title = "with B7, D7")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
dev.off()

write_csv(x = counts[counts$allele.species.consistency=="Consistent with genotype",]%>%
            mutate(allele = substr(x = allele,start = 1,stop = 4))%>%
            select(-allele.species.consistency)%>%
            group_by(lib,species,allele,treatment,tissue,dmel_gene,SC_sig_DE.cat)%>%
            summarise(counts= sum(counts)),file = "Starting.expression.data.for.edgeR.csv")

