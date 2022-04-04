library(tidyverse)
library(markdown)
library(ggtext)
library(lme4)

remove(list = ls(all.names = T))

miseq <- read.csv("miseq.csv")

pdf("miseq.allele.proportions_final.pdf",width = 6,height = 5)
ggplot(data = miseq)+
  geom_boxplot(aes(x = treatment,y = nuc.cover.prop,colour = allele),width = 0.5,outlier.shape = NA)+
  geom_jitter(aes(x = treatment,y = nuc.cover.prop,colour = allele),position = position_dodge2(0.4),size=0.5)+
  geom_hline(yintercept = 0.5,linetype=3,colour="grey")+
  # annotate(geom = "text",x = 4,y = 0.5,label="no~italic(cis)",parse=T,size=2.5)+
  # expand_limits(x=4.5)+
  facet_wrap(~seqnames,ncol = 5)+
  labs(x="",y="Proportions of transcripts")+
  # guides(colour=guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 0.9),
        # plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = c(0.85,0.1),
        legend.title = element_blank(),
        legend.text = element_markdown())
dev.off()



##########################################################################
#### Binomial regression  ####
##########################################################################
library(lme4)
remove(list = ls(all.names = T))
hyb.miseq <- read.csv("miseq.hyb.proportions.1.SNP.for.each.gene.csv")
hyb.miseq[is.na(hyb.miseq)]<-0
hyb.miseq.1 <- data.frame(gene=character(),lib=character(),treatment=character(),count.sim=numeric())
for (i in 1:length(hyb.miseq[,1])) {
  hyb.miseq.1<-bind_rows(hyb.miseq.1,
                         data.frame(gene=rep(hyb.miseq$seqnames[i],hyb.miseq$count.sim[i]),
                                    lib=rep(hyb.miseq$library[i],hyb.miseq$count.sim[i]),
                                    treatment=rep(hyb.miseq$treatment[i],hyb.miseq$count.sim[i]),
                                    count.sim=rep(1,hyb.miseq$count.sim[i])),
                         data.frame(gene=rep(hyb.miseq$seqnames[i],hyb.miseq$count.sec[i]),
                                    lib=rep(hyb.miseq$library[i],hyb.miseq$count.sec[i]),
                                    treatment=rep(hyb.miseq$treatment[i],hyb.miseq$count.sec[i]),
                                    count.sim=rep(0,hyb.miseq$count.sec[i])))
}

for (g in unique(hyb.miseq.1$gene)) {
  .GlobalEnv[[paste0("fit.",g)]]<-glmer(formula = count.sim~treatment+(1|lib),family = binomial,data = hyb.miseq.1%>%filter(treatment!="gDNA",gene==g))
  .GlobalEnv[[paste0("rslt.",g)]]<-summary(.GlobalEnv[[paste0("fit.",g)]])
}
result<-data.frame(gene=character(),
                   logit.ctrl=numeric(),logit.ctrl.CI.lower=numeric(),logit.ctrl.CI.upper=numeric(),
                   logit.wasp=numeric(),logit.wasp.CI.lower=numeric(),logit.wasp.CI.upper=numeric(),
                   p.intercept=numeric(),p.trt=numeric())
for (g in unique(hyb.miseq$seqnames)) {
  result<-result%>%add_row(
    gene = g,
    logit.ctrl = .GlobalEnv[[paste0("rslt.",g)]]$coefficients[1],
    logit.ctrl.CI.lower=.GlobalEnv[[paste0("rslt.",g)]]$coefficients[1]-.GlobalEnv[[paste0("rslt.",g)]]$coefficients[3]*1.96,
    logit.ctrl.CI.upper=.GlobalEnv[[paste0("rslt.",g)]]$coefficients[1]+.GlobalEnv[[paste0("rslt.",g)]]$coefficients[3]*1.96,
    
    logit.wasp=.GlobalEnv[[paste0("rslt.",g)]]$coefficients[1]+.GlobalEnv[[paste0("rslt.",g)]]$coefficients[2],
    logit.wasp.CI.lower=.GlobalEnv[[paste0("rslt.",g)]]$coefficients[1]+.GlobalEnv[[paste0("rslt.",g)]]$coefficients[2]-.GlobalEnv[[paste0("rslt.",g)]]$coefficients[4]*1.96,
    logit.wasp.CI.upper=.GlobalEnv[[paste0("rslt.",g)]]$coefficients[1]+.GlobalEnv[[paste0("rslt.",g)]]$coefficients[2]+.GlobalEnv[[paste0("rslt.",g)]]$coefficients[4]*1.96,
    
    p.intercept=.GlobalEnv[[paste0("rslt.",g)]]$coefficients[7],
    p.trt=.GlobalEnv[[paste0("rslt.",g)]]$coefficients[8]
  )
}
result.plot<-bind_rows(result%>%select(gene,logit.ctrl,logit.ctrl.CI.lower,logit.ctrl.CI.upper)%>%`colnames<-`(c("gene","logit","logit.CI.lower","logit.CI.upper"))%>%mutate(treatment='Control'),
                       result%>%select(gene,logit.wasp,logit.wasp.CI.lower,logit.wasp.CI.upper)%>%`colnames<-`(c("gene","logit","logit.CI.lower","logit.CI.upper"))%>%mutate(treatment='Wasp'))%>%mutate(gene=factor(x = gene,levels = c("BomBc1","BomBc2","BomBc3","BomS1",
                                                                                                                                                                                                                                         "Mtk","CecA1","CecC","edin",
                                                                                                                                                                                                                                         "TotB","TotX","CG16772","CG8850","Indy")))
pdf("miseq.hyb.allele.logit.pdf",width = 6,height = 3)
ggplot(data = result.plot,mapping = aes(x = gene,y = logit,colour=treatment))+
  scale_color_manual(values = list("Control"="grey","Wasp"="black"))+
  geom_errorbar(mapping = aes(ymin=logit.CI.lower,ymax=logit.CI.upper),width=.3,position = position_dodge(0.5))+
  geom_point(position = position_dodge(0.5))+
  geom_hline(yintercept = 0,linetype=2,colour='pink')+
  expand_limits(x=14.5)+
  annotate(geom = "text",x = 13.8,y = 0.2,label="no~italic(cis)",parse=T,size=4)+
  # scale_y_continuous(breaks = c(-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5))+
  # scale_y_continuous(breaks = c(-3.5,-2.5,-1.5,-0.5,0,0.5,1.5,2.5,3.5))+
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3))+
  ylab("Logit( *D. simulans* allele )")+
  guides(colour=guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 0.9),
        legend.title = element_blank(),
        # plot.margin = unit(c(1.3,0.5,0.5,0.5),"cm"),
        legend.position = c(0.85,0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown())
dev.off()


