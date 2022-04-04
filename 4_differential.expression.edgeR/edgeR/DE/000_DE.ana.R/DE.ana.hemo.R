library(edgeR)

library(tidyverse)
library(ggrepel)
library(fdrtool)

rm(list = ls(all.names = TRUE))

y.counts <- read.csv("edgeR.DE.data.csv",header = T,row.names = 1)
y.counts <- y.counts[,grepl(pattern = "hemo",x = colnames(y.counts))]
colnames(y.counts)<-gsub(pattern = ".hemo",replacement = "",x = colnames(y.counts))
# import genes:
genes <- read.table("ASDE.analysed.genes.for.haemocytes.txt",col.names = "gene")





###################################################################
####################################### creating DGEList ##########
###################################################################
y <- DGEList(counts = y.counts,group = substr(x=names(y.counts),start = 5,stop = 17))
keep <- row.names(y$counts)%in%c(genes$gene)
table(keep)
# keep
# FALSE  TRUE 
# 2049  5275 
y<-y[keep, , keep.lib.sizes=FALSE]
y<-calcNormFactors(y)
### MD plot: column=sample(i)
plotMD(cpm(y, log=TRUE), column=2)
abline(h=0, col="red", lty=2, lwd=2)
write.csv(x = cpm(y),file = "within.cpm.hemo.csv")
###################################################################
################# design ##########################################
###################################################################
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
### Estimating the dispersion:
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
# [1] 0.1417077
plotBCV(y)
### fit quasi-likelihood model
fit<-glmQLFit(y,design = design,robust=TRUE)
plotQLDisp(fit)


#############################
#####  ################
#############################
### set contrast:
con.DE <- makeContrasts(
  con.DE.sim = SIM.wasp-SIM.ctrl,
  con.DE.sec = SEC.wasp-SEC.ctrl,
  con.DE.hyb = HYB.wasp-HYB.ctrl,
  
  con.DE.h.sim.c = HYB.ctrl - SIM.ctrl,
  con.DE.h.sim.w = HYB.wasp - SIM.wasp,
  con.DE.h.sim.r = (HYB.wasp-HYB.ctrl) - (SIM.wasp-SIM.ctrl),
  con.DE.h.sec.c = HYB.ctrl - SEC.ctrl,
  con.DE.h.sec.w = HYB.wasp - SEC.wasp,
  con.DE.h.sec.r = (HYB.wasp-HYB.ctrl) - (SEC.wasp-SEC.ctrl),
  con.DE.sim.sec.c = SIM.ctrl - SEC.ctrl,
  con.DE.sim.sec.w = SIM.wasp - SEC.wasp,
  con.DE.sim.sec.r = (SIM.wasp-SIM.ctrl) - (SEC.wasp-SEC.ctrl),
  
  levels = design
)
### QL F-test:
for (con in colnames(con.DE)) {
  .GlobalEnv[[paste0('QLFT_DE.',substr(x = con,start = 8,stop = nchar(con)))]] <- glmQLFTest(fit,contrast = con.DE[,con])
  .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]] <- topTags(object = .GlobalEnv[[paste0('QLFT_DE.',substr(x = con,start = 8,stop = nchar(con)))]],n = Inf)$table
  .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]<-.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]%>%
    mutate(dmel_gene = rownames(.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]))
  rownames(.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]])<-c()
}


### arrnage results:

result.within.DE <- bind_rows(df.result_DE.sim%>%mutate(species='*D. simulans*'),
                              df.result_DE.sec%>%mutate(species='*D. sechellia*'),
                              df.result_DE.hyb%>%mutate(species='Hybrid'))%>%mutate(tissue='Haemocytes')
write.csv(x = result.within.DE, file = "DE.result.hemo.response.within.SIM.SEC.HYB.csv",row.names = F)

result.h.sim <- bind_rows(df.result_DE.h.sim.c%>%mutate(treatment='Control'),
                          df.result_DE.h.sim.w%>%mutate(treatment='Wasp'),
                          df.result_DE.h.sim.r%>%mutate(treatment='Response'))%>%mutate(tissue='Haemocytes')
result.h.sec <- bind_rows(df.result_DE.h.sec.c%>%mutate(treatment='Control'),
                          df.result_DE.h.sec.w%>%mutate(treatment='Wasp'),
                          df.result_DE.h.sec.r%>%mutate(treatment='Response'))%>%mutate(tissue='Haemocytes')
result.between.DE <- merge.data.frame(x = result.h.sim,y = result.h.sec,by = c('dmel_gene','tissue','treatment'),suffixes = c('.hyb_sim','.hyb_sec'))
write.csv(x = result.between.DE, file = "DE.result.hemo.response.betweenHSim.betweenHsec.csv",row.names = F)

result.sim.sec <- bind_rows(df.result_DE.sim.sec.c%>%mutate(treatment='Control'),
                            df.result_DE.sim.sec.w%>%mutate(treatment='Wasp'),
                            df.result_DE.sim.sec.r%>%mutate(treatment='Response'))%>%mutate(tissue='Haemocytes')
write.csv(x = result.sim.sec, file = "DE.result.hemo.response.betweenSIMSEC.csv",row.names = F)












