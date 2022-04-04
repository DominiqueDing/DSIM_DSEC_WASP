library(ggplot2)
library(tidyverse)
library(ggsignif)
library(Hmisc)
library(ggtext)
library(ggpubr)
library(lemon)


#########################################
#### statistic analyses   ###############
#########################################
remove(list = ls())
load("CT.weighted.RData")
data.0 <- CT.weighted%>%select(-x)%>%drop_na()%>%filter(!gene%in%c("rpl","rpl_new"))%>%
  mutate(allele=factor(allele,levels = c("*D. sechellia*","*D. simulans*")),
         treatment=factor(treatment,levels = c("Control","Wasp")),
         species=factor(species,levels = c("Hybrid","Parents")))

############## For plotting:
data.0$x <- paste0(data.0$species,'_',data.0$allele,'_',data.0$treatment)
for (g in unique(data.0$gene)) {
  .GlobalEnv[[paste0("fit.plot.",g)]] <- lm(ct ~ 0 + x,data = data.0%>%filter(gene==g))
  .GlobalEnv[[paste0("sum.plot.",g)]] <- summary(.GlobalEnv[[paste0("fit.plot.",g)]])
}

g <- unique(data.0$gene)
result.plot <- as.data.frame(sum.plot.BomBc1$coefficients)%>%mutate(x=rownames(.),gene=g[1])%>%`rownames<-`(.,NULL)
for (i in 2:length(g)) {
  result.plot <- bind_rows(
    result.plot,
    as.data.frame(.GlobalEnv[[paste0("sum.plot.",g[i])]]$coefficients)%>%mutate(x=rownames(.),gene=g[i])%>%`rownames<-`(.,NULL)
    )
}
result.plot<-result.plot%>%mutate(species=sub(pattern = "x",replacement = "",x = sub(pattern = "_.*",replacement ="",x = x)),
                                  allele=sub(pattern = "_.*",replacement = "",x = sub(pattern = ".*?_",replacement = "",x = x)),
                                  treatment=sub(pattern = ".*_",replacement = "",x = sub(pattern = ".*_",replacement ="",x = x)),
                                  allele=factor(allele,levels = c("*D. simulans*","*D. sechellia*")))
result.plot$treatment[result.plot$species=="Parents"]<-paste0(" ",result.plot$treatment[result.plot$species=="Parents"])
result.plot$group<-0
result.plot$group[result.plot$species=="Parents"&result.plot$treatment==" Control"&result.plot$allele=="*D. simulans*"]<-1
result.plot$group[result.plot$species=="Parents"&result.plot$treatment==" Wasp"&result.plot$allele=="*D. simulans*"]<-1
result.plot$group[result.plot$species=="Parents"&result.plot$treatment==" Control"&result.plot$allele=="*D. sechellia*"]<-2
result.plot$group[result.plot$species=="Parents"&result.plot$treatment==" Wasp"&result.plot$allele=="*D. sechellia*"]<-2
result.plot$group[result.plot$species=="Hybrid"&result.plot$treatment=="Control"&result.plot$allele=="*D. simulans*"]<-3
result.plot$group[result.plot$species=="Hybrid"&result.plot$treatment=="Wasp"&result.plot$allele=="*D. simulans*"]<-3
result.plot$group[result.plot$species=="Hybrid"&result.plot$treatment=="Control"&result.plot$allele=="*D. sechellia*"]<-4
result.plot$group[result.plot$species=="Hybrid"&result.plot$treatment=="Wasp"&result.plot$allele=="*D. sechellia*"]<-4
result.plot <- merge.data.frame(x = result.plot,y = result.plot%>%group_by(gene)%>%summarise(ytop = max(Estimate+1.96*`Std. Error`)),by = c("gene"),all = T)
result.plot <- merge.data.frame(x = result.plot,y = result.plot%>%group_by(gene)%>%summarise(ybot = min(Estimate-1.96*`Std. Error`)),by = c("gene"),all = T)
result.plot<-result.plot%>%mutate(gene=paste0("*",gene,"*"))
result.plot$gene[grepl("CG",result.plot$gene)]<-gsub(pattern = "\\*",replacement = '',result.plot$gene[grepl("CG",result.plot$gene)])

pdf("fatbody_qpcr_final.pdf",width = 7,height = 7.5)
ggplot(data = result.plot,
       mapping = aes(x = treatment,
                     y = Estimate,
                     colour=allele))+
  geom_point(position = position_dodge(0.3))+
  geom_errorbar(mapping = aes(ymin=Estimate-1.96*`Std. Error`,ymax=Estimate+1.96*`Std. Error`,colour=allele),width=0.3,position = position_dodge(0.3))+
  geom_line(aes(group=group))+
  geom_hline(aes(yintercept = ytop+(ytop-ybot)/10),colour="dark grey")+
  geom_hline(aes(yintercept = ytop+(ytop-ybot)/3),colour="white")+
  geom_vline(xintercept = 2.5,colour="dark grey")+
  geom_text(aes(x = 1.5,y = ytop+(ytop-ybot)/4,label=ifelse(test = species=="Parents",yes = species,no = "")),colour=adjustcolor( "black", alpha.f = 0.3),size=3)+
  geom_text(aes(x = 3.5,y = ytop+(ytop-ybot)/4,label=ifelse(test = species=="Hybrid",yes = species,no = "")),colour=adjustcolor( "black", alpha.f = 0.3),size=3)+
  facet_wrap(facets = ~gene,ncol = 4,scales = 'free_y')+
  ylab(label = "Log<sub>2</sub>( Relative gene expression )")+
  guides(colour=guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(legend.text = element_markdown(size=10),
        legend.title = element_blank(),
        legend.position = c(0.5,0.1),
        panel.grid = element_blank(),
        # panel.spacing = unit(0,'lines'),
        strip.background = element_blank(),
        strip.text = element_markdown(size = 10),
        axis.title.y = element_markdown(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 0.9))
dev.off()
























#########################################
#### Estimate trans effects   ###########
#########################################

library(lmtest)
#### is trans significant?
for (g in unique(data.0$gene)) {
  .GlobalEnv[[paste0("full.",g)]] <- lm(ct ~ allele*species*treatment,data = data.0%>%filter(gene==g))
  .GlobalEnv[[paste0("sum.full.",g)]] <- summary(.GlobalEnv[[paste0("full.",g)]])
  
  .GlobalEnv[[paste0("sans.trans.",g)]] <- lm(ct ~ allele + treatment +species +allele:treatment,data = data.0%>%filter(gene==g))
  .GlobalEnv[[paste0("sum.sans.trans.",g)]] <- summary(.GlobalEnv[[paste0("sans.trans.",g)]])
  
  .GlobalEnv[[paste0("lrt.1.",g)]] <- lrtest(.GlobalEnv[[paste0("full.",g)]],.GlobalEnv[[paste0("sans.trans.",g)]])
  .GlobalEnv[[paste0("sum.lrt.1.",g)]] <- summary(.GlobalEnv[[paste0("lrt.1.",g)]])
}
result.trans.sig <- data.frame(gene=unique(data.0$gene))
result.trans.sig$P <- 0
for (i in 1:length(result.trans.sig$gene)) {
  result.trans.sig$P[i] <- .GlobalEnv[[paste0("lrt.1.",result.trans.sig$gene[i])]]$`Pr(>Chisq)`[2]
}




library(alr4)
#### Trans-effect under control condition:
g <- unique(data.0$gene)
result.trans.plot.control <- as.data.frame(sum.full.BomBc1$coefficients)%>%mutate(x=rownames(.),gene=g[1])%>%`rownames<-`(.,NULL)%>%filter(x=="allele*D. simulans*:speciesParents")
for (i in 2:length(g)) {
  result.trans.plot.control <- bind_rows(
    result.trans.plot.control,
    as.data.frame(.GlobalEnv[[paste0("sum.full.",g[i])]]$coefficients)%>%mutate(x=rownames(.),gene=g[i])%>%`rownames<-`(.,NULL)%>%filter(x=="allele*D. simulans*:speciesParents")
  )
}
result.trans.plot.control$treatment <- "Control"
#### Trans-effect under infected condition:
for (g in unique(data.0$gene)) {
  .GlobalEnv[[paste0("full.0.1.",g)]] <- lm(ct ~ allele*species*treatment,data = data.0%>%filter(gene==g)%>%mutate(treatment=factor(treatment,levels=c("Wasp","Control"))))
  .GlobalEnv[[paste0("sum.full.0.1.",g)]] <- summary(.GlobalEnv[[paste0("full.0.1.",g)]])
}
g <- unique(data.0$gene)
result.trans.plot.wasp <- as.data.frame(sum.full.0.1.BomBc1$coefficients)%>%mutate(x=rownames(.),gene=g[1])%>%`rownames<-`(.,NULL)%>%filter(x=="allele*D. simulans*:speciesParents")
for (i in 2:length(g)) {
  result.trans.plot.wasp <- bind_rows(
    result.trans.plot.wasp,
    as.data.frame(.GlobalEnv[[paste0("sum.full.0.1.",g[i])]]$coefficients)%>%mutate(x=rownames(.),gene=g[i])%>%`rownames<-`(.,NULL)%>%filter(x=="allele*D. simulans*:speciesParents")
  )
}
result.trans.plot.wasp$treatment <- "Wasp"




result.trans.plot <- bind_rows(result.trans.plot.control,result.trans.plot.wasp)

pdf("qpcr.trans.plot.pdf",width = 6,height = 3)
ggplot(data = result.trans.plot,
       mapping = aes(x = gene,y = Estimate,colour=treatment))+
  scale_color_manual(values = list("Control"="grey","Wasp"="black"))+
  geom_errorbar(mapping = aes(ymin=Estimate-1.96*`Std. Error`,ymax=Estimate+1.96*`Std. Error`),width=.3,position = position_dodge(0.5))+
  geom_point(position = position_dodge(0.5))+
  geom_hline(yintercept = 0,linetype=2,colour='pink')+
  expand_limits(x=15.3)+
  annotate(geom = "text",x = 14.4,y = 0.32,label="no~italic(trans)",parse=T,size=4)+
  # scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3))+
  ylab("*Trans*-effect (Log<sub>2</sub> relative expression)")+
  guides(colour=guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 0.9),
        legend.title = element_blank(),
        plot.margin = unit(c(1.3,0.5,0.5,0.5),"cm"),
        legend.position = c(0.85,1.1),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=10))
dev.off()




#### Estimates of interaction:
result.trans.interaction <- data.frame(gene=unique(data.0$gene))
result.trans.interaction$P <- 0
result.trans.interaction$est <- 0
for (i in 1:length(result.trans.interaction$gene)) {
  result.trans.interaction$est[i] <- .GlobalEnv[[paste0("sum.full.",result.trans.interaction$gene[i])]]$coefficients[8]
  result.trans.interaction$P[i] <- .GlobalEnv[[paste0("sum.full.",result.trans.interaction$gene[i])]]$coefficients[32]
}


