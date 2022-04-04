library(ggplot2)
library(tidyverse)
library(ggsignif)
library(Hmisc)
library(ggtext)
library(ggpubr)
library(lemon)

load(file = "mel_assay.RData")


### calculate melinasation point estimates, upper.CI and lower.CI 95% CI for each species
testFunc <- function(a,b) binconf(x = a,n = a+b,method = "wilson")
apply(mel.adults[,c('mel','non.mel')], 1, function(x) testFunc(x[1],x[2]))

mel.adults <- bind_cols(mel.adults,
                        as.data.frame(t(apply(mel.adults[,c('mel','non.mel')], 1, function(x) testFunc(x[1],x[2])))))
names(mel.adults) <- c(names(mel.adults)[1:3],'mean.mel.rate','lower.CI','upper.CI')
mel.adults[1,4:6] <- 0
mel.adults<-mel.adults%>%mutate(species = factor(species,
                                                 levels = c("*D. simulans*","Hybrid","*D. sechellia*")))%>%arrange(species)

mel.larva <- bind_cols(mel.larva,
                        as.data.frame(t(apply(mel.larva[,c('mel','non.mel')], 1, function(x) testFunc(x[1],x[2])))))
names(mel.larva) <- c(names(mel.larva)[1:3],'mean.mel.rate','lower.CI','upper.CI')
mel.larva[1,4:6] <- 0
mel.larva<-mel.larva%>%arrange()%>%mutate(species = factor(species,
                                                 levels = c("*D. simulans*","Hybrid","*D. sechellia*")))%>%arrange(species)


### Binomial regression
model.adults<- glm(formula = cbind(mel.adults$mel[1:2],mel.adults$non.mel[1:2])~mel.adults$species[1:2],
                      family = "binomial")
summary(model.adults)

model.larva<- glm(formula = cbind(mel.larva$mel[1:2],mel.larva$non.mel[1:2])~mel.larva$species[1:2],
                   family = "binomial")
summary(model.larva)


### plot
pdf("../figures/injec_assay_larv.pdf",width = 3,height = 5)
ggplot(data = mel.larva,
       mapping = aes(x = species,y = mean.mel.rate))+
  geom_col(fill="dark grey",width = 0.5)+
  geom_errorbar(mapping = aes(ymin=lower.CI,ymax=upper.CI),
                width=0.2)+
  annotate(geom = "text",label = "N = 122", x = 1,y = 0.16)+
  annotate(geom = "text",label = "N = 102", x = 2,y = 0.98)+
  annotate(geom = "text",label = "", x = 2,y = 1.05)+
  annotate(geom = "text",label = "N = 72", x = 3,y = 0.05)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),limits = c(0,1.05),expand = c(0.05,0))+
  labs(y = "Rate Oil Droplets Encapsulated (24hpi)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y= element_text(colour="black",size=14),
        axis.text.x = element_markdown(colour="black",size=14,angle = 45,vjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_blank())+
  coord_capped_cart(left = capped_vertical(capped = "top"))
# +
#   geom_signif(comparisons = list(c("*D. simulans*","Hybrid"),
#                                  c("*D. sechellia*","Hybrid"),
#                                  c("*D. sechellia*","*D. simulans*")),
#               y_position = c(1., 1.05, 1.12),
#               annotations = c("***","***","***"))
dev.off()
pdf("../figures/injec_assay_fly.pdf",width = 3,height = 5)
ggplot(data = mel.adults,
       mapping = aes(x = species,y = mean.mel.rate))+
  geom_col(fill="dark grey",width = 0.5)+
  geom_errorbar(mapping = aes(ymin=lower.CI,ymax=upper.CI),
                width=0.2)+
  annotate(geom = "text",label = "N = 88", x = 1,y = 0.8)+
  annotate(geom = "text",label = "N = 83", x = 2,y = 1.05)+
  annotate(geom = "text",label = "N = 63", x = 3,y = 0.05)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),limits = c(0,1.05),expand = c(0.05,0))+
  labs(y = "Rate Oil Droplets Encapsulated (adults)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y= element_text(colour="black",size=14),
        axis.text.x = element_markdown(colour="black",size=14,angle = 45,vjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_blank())+
  coord_capped_cart(left = capped_vertical(capped = "top"))
# +
#   geom_signif(comparisons = list(c("*D. simulans*","Hybrid"),
#                                  c("*D. sechellia*","Hybrid"),
#                                  c("*D. sechellia*","*D. simulans*")),
#               y_position = c(1.05, 1.1, 1.17),
#               annotations = c("**","***","***"))

dev.off()


