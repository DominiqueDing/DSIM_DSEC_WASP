# library("edgeR", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(edgeR)
library(tidyverse)
library(readxl)

remove(list = ls(all.names = T))

# load data prepared for edgeR analyses: (-mismap>0.1 genes; dsim,dsec true to allele reads)
counts <- read_csv("Starting.expression.data.for.edgeR.csv",col_names = T)

# prepare dataset suitable for edgeR analyses:
counts.y<-counts%>%select(-SC_sig_DE.cat,-allele)
# combine allele-specific reads in hybrid libraries:
counts.y <- counts.y%>%group_by(lib,species,treatment,tissue,dmel_gene)%>%summarise(counts=sum(counts))


counts.y$species<-sub(pattern = "\\*D. sechellia\\*",replacement = "SEC",x = counts.y$species)
counts.y$species<-sub(pattern = "\\*D. simulans\\*",replacement = "SIM",x = counts.y$species)
counts.y$species<-sub(pattern = "Hybrid",replacement = "HYB",x = counts.y$species)

counts.y$treatment<-sub(pattern = "Control",replacement = "ctrl",x = counts.y$treatment)
counts.y$treatment<-sub(pattern = "Wasp",replacement = "wasp",x = counts.y$treatment)

counts.y$tissue<-sub(pattern = "Haemocytes",replacement = "hemo",x = counts.y$tissue)
counts.y$tissue<-sub(pattern = "Fat body",replacement = "fabo",x = counts.y$tissue)

counts.y$counts.ID<-paste(counts.y$lib,counts.y$species,counts.y$treatment,counts.y$tissue,sep = '.')
counts.y$group.ID<-paste(counts.y$species,counts.y$treatment,counts.y$tissue,sep = '.')


counts.y.lst <- list()
dmel_gene<-unique(counts.y$dmel_gene)
counts.ID<-unique(counts.y$counts.ID)
for (i in 1:length(counts.ID)) {
  counts.y.lst[[i]] <- counts.y[counts.y$counts.ID==counts.ID[i],c("dmel_gene","counts")]
  names(counts.y.lst[[i]])<-c("dmel_gene",counts.ID[i])
}
y.counts<-counts.y.lst[[1]]
for (i in 2:length(counts.y.lst)) {
  y.counts<-merge(x = y.counts,y = counts.y.lst[[i]],by = "dmel_gene")
}
remove(counts.y.lst)
row.names(y.counts)<-y.counts$dmel_gene
y.counts<-y.counts[,!names(y.counts)%in%c("dmel_gene")]

### proof editting y.counts names ###
for (bar in c("A","B","C","D")) {
  for (i in 1:9) {
    names(y.counts)<-gsub(pattern = paste(bar,i,"\\.",sep = ""),replacement = paste(bar,"0",i,"\\.",sep = ""),x = names(y.counts))
  }
}

write.csv(x = y.counts,file = "edgeR.DE.data.csv")
