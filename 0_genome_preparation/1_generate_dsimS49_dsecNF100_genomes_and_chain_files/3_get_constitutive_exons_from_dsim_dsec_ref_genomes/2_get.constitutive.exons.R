
library(tidyverse)

dsecREF.gene.mRNA <- read.table("dsecREF.mRNA.gene.tab",col.names = c("mRNA","gene"),fill = T)
dsecREF.gene.mRNA$mRNA.per.gene<-1
dsecREF.counts.mRNA.per.gene <- dsecREF.gene.mRNA%>%group_by(gene)%>%summarise(mRNA.per.gene=sum(mRNA.per.gene))

dsecREF.exon.mRNA <- read.table("dsecREF.exon.mRNA.tab",col.names = c("exon","mRNA"),fill = T)
dsecREF.exon.mRNA$mRNA.per.exon <-1
dsecREF.exon.mRNA<-merge(x = dsecREF.exon.mRNA,y = dsecREF.gene.mRNA,by = "mRNA")
dsecREF.counts.mRNA.per.exon <- dsecREF.exon.mRNA%>%group_by(exon,gene)%>%summarise(mRNA.per.exon=sum(mRNA.per.exon))
dsecREF.exon.gene<-merge(x = dsecREF.counts.mRNA.per.exon,y = dsecREF.counts.mRNA.per.gene,by = "gene",all = T)
dsecREF.constitutive.exon<-dsecREF.exon.gene[dsecREF.exon.gene$mRNA.per.exon==dsecREF.exon.gene$mRNA.per.gene,]
write(x=as.character(dsecREF.constitutive.exon$exon),file = "dsecREF.constitutive.exons.txt")










dsimREF.mRNA.gene <- read.table("dsimREF.mRNA.gene.tab",col.names = c("mRNA","gene"),fill = T)
dsimREF.mRNA.gene$mRNA.per.gene<-1
dsimREF.counts.mRNA.per.gene <- dsimREF.mRNA.gene%>%group_by(gene)%>%summarise(mRNA.per.gene=sum(mRNA.per.gene))
# sanity check
#dsimREF.mRNA.gene$gene.per.mRNA<-1
#dsimREF.counts.gene.per.mRNA <- dsimREF.mRNA.gene%>%group_by(mRNA)%>%summarise(gene.per.mRNA=sum(gene.per.mRNA))

dsimREF.exon.mRNA <- read.table("dsimREF.exon.mRNA.tab",col.names = c("exon",paste0("mRNA",1:44)),na.strings=c("","NA"),fill = T)
dsimREF.exon.mRNA.lst<-list()
for (i in 2:45) {
  dsimREF.exon.mRNA.lst[[i-1]]<-`colnames<-`(dsimREF.exon.mRNA[,c(1,i)],c("exon","mRNA"))%>%drop_na(mRNA)
}
dsimREF.exon.mRNA<-dsimREF.exon.mRNA.lst[[1]]
for (i in 2:length(dsimREF.exon.mRNA.lst)) {
  dsimREF.exon.mRNA<-bind_rows(dsimREF.exon.mRNA,dsimREF.exon.mRNA.lst[[i]])
}
dsimREF.exon.mRNA$mRNA.per.exon<-1
dsimREF.exon.mRNA<-merge(x = dsimREF.exon.mRNA,y=dsimREF.mRNA.gene,by = "mRNA")
dsimREF.counts.mRNA.per.exon<-dsimREF.exon.mRNA%>%group_by(exon,gene)%>%summarise(mRNA.per.exon=sum(mRNA.per.exon))
dsimREF.exon.gene<-merge(x = dsimREF.counts.mRNA.per.exon,y = dsimREF.counts.mRNA.per.gene,by = "gene",all = T)
dsimREF.constitutive.exon<-dsimREF.exon.gene[dsimREF.exon.gene$mRNA.per.exon==dsimREF.exon.gene$mRNA.per.gene,]
write(x = as.character(dsimREF.constitutive.exon$exon),file = "dsimREF.constitutive.exons.txt")
