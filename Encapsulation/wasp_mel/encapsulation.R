library(tidyverse)
library(readxl)
library(Hmisc)
library(ggtext)
library(lemon)


remove(list = ls())

wasp.assay <- read_excel("larvae.xlsx",sheet = "Sheet1")
wasp.assay.plot <- wasp.assay%>%group_by(species)%>%summarise(cap=sum(larvae_w_Capsules), wasp.larvae=sum(larvae_w_wasp_larvae))%>%mutate(N=cap+wasp.larvae)
### calculate melinasation point estimates, upper.CI and lower.CI 95% CI for each species
testFunc <- function(a,b) binconf(x = a,n = a+b,method = "wilson")
wasp.assay.plot <- bind_cols(wasp.assay.plot,
                        as.data.frame(t(apply(wasp.assay.plot[,c('cap','wasp.larvae')], 1, function(x) testFunc(x[1],x[2])))))
names(wasp.assay.plot) <- c(names(wasp.assay.plot)[1:4],'mean.cap.rate','lower.CI','upper.CI')
wasp.assay.plot[1,5:7] <- 0
# wasp.assay.plot[3,] <- rep(NA,7) 
# wasp.assay.plot[3,1] <- 'x'

pdf("../figures/wasp_assay.pdf",width = 2,height = 5)
ggplot(data = wasp.assay.plot%>%mutate(species=factor(species,levels = c("*D. simulans*","*D. sechellia*",'x'))),
       mapping = aes(x = species,y = mean.cap.rate))+
  geom_col(fill="dark grey",width=0.5)+
  geom_errorbar(mapping = aes(ymin=lower.CI,ymax=upper.CI),width=0.2)+
  annotate(geom = "text",label = "N = 92", x = 1,y = 0.3)+
  annotate(geom = "text",label = "N = 22", x = 2,y = 0.05)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),limits = c(0,1.05),expand = c(0.05,0))+coord_capped_cart(left = capped_vertical(capped = "top"))+
  # scale_y_continuous(breaks = c(0,0.15,0.25,0.5,0.75,1),expand = c(0.05,0))+
  scale_x_discrete(breaks = c("*D. simulans*","*D. sechellia*"),labels=c("*D. simulans*","*D. sechellia*"))+
  labs(y = "Rate Wasp Eggs Encapsulated (96hpi)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y= element_text(colour="black",size=14),
        axis.text.x = element_markdown(colour="black",size=14,angle = 45,vjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_blank())
dev.off()
