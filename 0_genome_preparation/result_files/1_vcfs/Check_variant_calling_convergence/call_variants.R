
library(tidyverse)

vcf.counts.dsim49<-read.table("dsim49.vcf.counts.txt")
vcf.counts.dsecNF100<-read.table("dsecNF100.vcf.counts.txt")
##### plot convergence
pdf("dsim49.BQSR.pdf")
ggplot(data = vcf.counts.dsecNF100,aes(x= c(1:length(V1)),y = V1))+
  geom_line()+
  labs(x="Round of BQSR",y="Number of variants",title = "Calling variants of DsecNF100 on Dsec_r1.3")
dev.off()
pdf("dsecNF100.BQSR.pdf")
ggplot(data = vcf.counts.dsim49,aes(x= c(1:length(V1)),y = V1))+
  geom_line()+
  labs(x="Round of BQSR",y="Number of variants",title = "Calling variants of Dsim49 on Dsim_r2.01")
dev.off()

GQ.dsim49 <- read.table("dsim49.vcf.counts/dsim49.round.all.BQSR.GQ.txt",col.names = c("GQ","round"),blank.lines.skip = T,na.strings = c("","NA"),sep = "\t")
GQ.dsecNF100 <- read.table("dsecNF100.vcf.counts/dsecNF100.round.all.BQSR.GQ.txt",col.names = c("GQ","round"),blank.lines.skip = T,na.strings = c("","NA"),sep = "\t")
GQ.dsim49<-na.omit(GQ.dsim49)
GQ.dsecNF100<-na.omit(GQ.dsecNF100)

### dsecNF100
# all histograms of GQ of every round:
pdf("dsecNF100.BQSR.GQ.histo.pdf")
ggplot(data = GQ.dsecNF100,aes(x=GQ))+
  geom_histogram(binwidth = 3)+
  facet_wrap(~round,scales = "free",ncol = 4)+
  labs(colour="",shape="",x="GQ score",y="",title = "GQ score histograms of each round of BQSR (DsecNF100)")
dev.off()
# histo compare 1st and 20th round
ggplot(data = as.data.frame(rbind(GQ.dsecNF100[GQ.dsecNF100$round==20,],GQ.dsecNF100[GQ.dsecNF100$round==1,])),
       aes(x=GQ,fill=as.factor(round)))+
  geom_histogram(binwidth=3,color="white",alpha=.3,position = "identity")
# frequency polygon compare 1st and 20th round
#ggplot(data = GQ.dsecNF100,aes(GQ,colour=as.factor(round)))+
#  geom_freqpoly(binwidth=3)
ggplot(data = rbind(GQ.dsecNF100[GQ.dsecNF100$round==1,],GQ.dsecNF100[GQ.dsecNF100$round==20,]),
       aes(GQ,colour=as.factor(round)))+
  geom_freqpoly(binwidth=3)

### dsim49
# all histograms of GQ of every round
pdf("dsim49.BQSR.GQ.histo.pdf")
ggplot(data = GQ.dsim49,aes(x=GQ))+
  geom_histogram(binwidth = 3)+
  facet_wrap(~round,scales = "free",ncol = 4)+
  labs(colour="",shape="",x="GQ score",y="",title = "GQ score histograms of each round of BQSR (Dsim49)")
dev.off()
# histogram compare every round
ggplot(data = GQ.dsim49,aes(x=GQ,fill=as.factor(round)))+
  geom_histogram(binwidth=3,colour="white",alpha=.3,position = "identity")
# histogram compare all round sans 1st
ggplot(data = GQ.dsim49[GQ.dsim49$round!=1,],aes(x=GQ,fill=as.factor(round)))+
  geom_histogram(binwidth=3,colour="white",alpha=.3,position = "identity")
# polygon compare all rounds
ggplot(data = GQ.dsim49,aes(GQ,colour=as.factor(round)))+
  geom_freqpoly(binwidth=3)
# polygon cmopare all rounds sans 1st
ggplot(data = GQ.dsim49[GQ.dsim49$round!=1,],aes(GQ,colour=as.factor(round)))+
  geom_freqpoly(binwidth=3)
# histogram compare 1st and 24th round:
ggplot(data = rbind(GQ.dsim49[GQ.dsim49$round==1,],GQ.dsim49[GQ.dsim49$round==24,]),
       aes(x=GQ,fill=as.factor(round)))+
  geom_histogram(binwidth=3,colour="white",alpha=.5,position="identity")
# polugon compare 1st and 24th round
ggplot(data = rbind(GQ.dsim49[GQ.dsim49$round==1,],GQ.dsim49[GQ.dsim49$round==20,]),
       aes(GQ,colour=as.factor(round)))+
  geom_freqpoly(binwidth=3)

