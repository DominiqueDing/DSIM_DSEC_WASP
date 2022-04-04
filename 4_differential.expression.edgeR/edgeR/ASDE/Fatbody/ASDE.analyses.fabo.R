# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")}
# BiocManager::install("edgeR")
library(edgeR)

library(tidyverse)
library(ggrepel)
library(fdrtool)

rm(list = ls(all.names = TRUE))






y.counts <- read.csv("edgeR.ASDE.data.csv",header = T,row.names = 1)
y.counts <- y.counts[,grepl(pattern = "fabo",x = colnames(y.counts))]
colnames(y.counts)<-gsub(pattern = "fabo.",replacement = "",x = colnames(y.counts))

###################################################################
####################################### creating DGEList ##########
###################################################################
y.fatbody <- DGEList(counts = y.counts,group = substr(x=names(y.counts),start = 5,stop = 17))
# remove(y.counts)
keep<-rowSums(cpm(y.fatbody)>3) >= 5
table(keep)
# keep
# FALSE  TRUE
#  2355  4969
y.fatbody<-y.fatbody[keep, , keep.lib.sizes=FALSE]
y.fatbody<-calcNormFactors(y.fatbody)
### MD plot: column=sample(i)
plotMD(cpm(y.fatbody, log=TRUE), column=2)
abline(h=0, col="red", lty=2, lwd=2)
write.table(x = rownames(y.fatbody$counts),file = "ASDE.analysed.genes.for.fatbody.txt",quote = F,row.names = F,col.names = F)

###################################################################
################# design ##########################################
###################################################################
design <- model.matrix(~ 0 + group, data = y.fatbody$samples)
colnames(design) <- levels(y.fatbody$samples$group)
### Estimating the dispersion:
y.fatbody <- estimateDisp(y.fatbody, design, robust = TRUE)
y.fatbody$common.dispersion
# [1] 0.1006463
plotBCV(y.fatbody)
### fit quasi-likelihood model
fit<-glmQLFit(y.fatbody,design = design,robust=TRUE)
plotQLDisp(fit)






#########################################
#####                               #####
#####    Controls (SIM/SEC/HYB)     #####
#####                               #####
#########################################
### set contrast:
con.ctrl <- makeContrasts(
  con.ctrl.total = SIM.ctrl.dsim - SEC.ctrl.dsec,
  con.ctrl.cis = HYB.ctrl.dsim - HYB.ctrl.dsec,
  con.ctrl.trans = (SIM.ctrl.dsim - SEC.ctrl.dsec) - (HYB.ctrl.dsim - HYB.ctrl.dsec),

  levels = design
)
### QL F-test:
for (con in colnames(con.ctrl)) {
  .GlobalEnv[[paste0('QLFT_ctrl.',substr(x = con,start = 10,stop = nchar(con)))]] <- glmQLFTest(fit,contrast = con.ctrl[,con])
  .GlobalEnv[[paste0('df.result_ctrl.',substr(x = con,start = 10,stop = nchar(con)))]] <- topTags(object = .GlobalEnv[[paste0('QLFT_ctrl.',substr(x = con,start = 10,stop = nchar(con)))]],n = Inf)$table
  .GlobalEnv[[paste0('df.result_ctrl.',substr(x = con,start = 10,stop = nchar(con)))]]<-.GlobalEnv[[paste0('df.result_ctrl.',substr(x = con,start = 10,stop = nchar(con)))]]%>%
    mutate(dmel_gene = rownames(.GlobalEnv[[paste0('df.result_ctrl.',substr(x = con,start = 10,stop = nchar(con)))]]))
}

### arrnage results:
ana.ctrl <- merge(x = df.result_ctrl.cis,y = df.result_ctrl.trans,by = "dmel_gene",suffixes = c(".cis",".trans"))
names(df.result_ctrl.total)[-6] <- paste0(names(df.result_ctrl.total)[-6],".total")
ana.ctrl <- merge(x = ana.ctrl,y = df.result_ctrl.total,by = "dmel_gene")

#########################################
#####                               #####
#####   Wasp (SIM/SEC/HYB)          #####
#####                               #####
#########################################
### set contrast:
con.wasp <- makeContrasts(
  con.wasp.total = SIM.wasp.dsim - SEC.wasp.dsec,
  con.wasp.cis = HYB.wasp.dsim - HYB.wasp.dsec,
  con.wasp.trans = (SIM.wasp.dsim - SEC.wasp.dsec) - (HYB.wasp.dsim - HYB.wasp.dsec),
  
  levels = design
)
### QL F-test:
for (con in colnames(con.wasp)) {
  .GlobalEnv[[paste0('QLFT_wasp.',substr(x = con,start = 10,stop = nchar(con)))]] <- glmQLFTest(fit,contrast = con.wasp[,con])
  .GlobalEnv[[paste0('df.result_wasp.',substr(x = con,start = 10,stop = nchar(con)))]] <- topTags(object = .GlobalEnv[[paste0('QLFT_wasp.',substr(x = con,start = 10,stop = nchar(con)))]],n = Inf)$table
  .GlobalEnv[[paste0('df.result_wasp.',substr(x = con,start = 10,stop = nchar(con)))]]<-.GlobalEnv[[paste0('df.result_wasp.',substr(x = con,start = 10,stop = nchar(con)))]]%>%
    mutate(dmel_gene = rownames(.GlobalEnv[[paste0('df.result_wasp.',substr(x = con,start = 10,stop = nchar(con)))]]))
}

### arrnage results:
ana.wasp <- merge(x = df.result_wasp.cis,y = df.result_wasp.trans,by = "dmel_gene",suffixes = c(".cis",".trans"))
names(df.result_wasp.total)[-6] <- paste0(names(df.result_wasp.total)[-6],".total")
ana.wasp <- merge(x = ana.wasp,y = df.result_wasp.total,by = "dmel_gene")

#############################################
#####                                   #####
#####  reac (SIM/SEC/HYB)               #####
#####                                   #####
#############################################
### set contrast:
con.reac <- makeContrasts(
  con.reac.total= (SIM.wasp.dsim-SIM.ctrl.dsim) - (SEC.wasp.dsec-SEC.ctrl.dsec),
  con.reac.cis = (HYB.wasp.dsim - HYB.ctrl.dsim) - (HYB.wasp.dsec - HYB.ctrl.dsec),
  con.reac.trans = ((SIM.wasp.dsim-SIM.ctrl.dsim) - (SEC.wasp.dsec-SEC.ctrl.dsec)) - ((HYB.wasp.dsim - HYB.ctrl.dsim) - (HYB.wasp.dsec - HYB.ctrl.dsec)),
  
  levels = design
)
### QL F-test:
for (con in colnames(con.reac)) {
  .GlobalEnv[[paste0('QLFT_reac.',substr(x = con,start = 10,stop = nchar(con)))]] <- glmQLFTest(fit,contrast = con.reac[,con])
  .GlobalEnv[[paste0('df.result_reac.',substr(x = con,start = 10,stop = nchar(con)))]] <- topTags(object = .GlobalEnv[[paste0('QLFT_reac.',substr(x = con,start = 10,stop = nchar(con)))]],n = Inf)$table
  .GlobalEnv[[paste0('df.result_reac.',substr(x = con,start = 10,stop = nchar(con)))]]<-.GlobalEnv[[paste0('df.result_reac.',substr(x = con,start = 10,stop = nchar(con)))]]%>%
    mutate(dmel_gene = rownames(.GlobalEnv[[paste0('df.result_reac.',substr(x = con,start = 10,stop = nchar(con)))]]))
}

### arrnage results:
ana.reac <- merge(x = df.result_reac.cis,y = df.result_reac.trans,by = "dmel_gene",suffixes = c(".cis",".trans"))
names(df.result_reac.total)[-6] <- paste0(names(df.result_reac.total)[-6],".total")
ana.reac <- merge(x = ana.reac,y = df.result_reac.total,by = "dmel_gene")


##############################
### category: ################
##############################
for (df in ls(pattern = "ana.")) {
  .GlobalEnv[[df]]$Category<-""
  .GlobalEnv[[df]]$Category.factor<-""
  
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<=0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<=0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Cis only"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Cis only"]<-1
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<=0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Ambiguous_cis.but.no.total"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_cis.but.no.total"]<-1.1
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<=0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(abs(.GlobalEnv[[df]]$logFC.cis-.GlobalEnv[[df]]$logFC)-
                                abs(.GlobalEnv[[df]]$logFC.trans-.GlobalEnv[[df]]$logFC)<=-1)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Ambiguous_total.but.no.cis.no.trans.yet.large.cisFC"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.cis.no.trans.yet.large.cisFC"]<-1.2 # large FC
  
  
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<=0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<=0.05)]<-"Trans only"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Trans only"]<-2
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<=0.05)]<-"Ambiguous_trans.but.no.total"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_trans.but.no.total"]<-2.1
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<=0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)
                            &(abs(.GlobalEnv[[df]]$logFC.trans-.GlobalEnv[[df]]$logFC.total)-
                                abs(.GlobalEnv[[df]]$logFC.cis-.GlobalEnv[[df]]$logFC.total)<=-1)]<-"Ambiguous_total.but.no.trans.no.cis.yet.large.transFC"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.trans.no.cis.yet.large.transFC"]<-2.2 # large FC
  
  
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<=0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<=0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<=0.05)
                            &(.GlobalEnv[[df]]$logFC.cis*.GlobalEnv[[df]]$logFC.trans<0)]<-"CisxTrans"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="CisxTrans"]<-3
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<=0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<=0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<=0.05)
                            &(.GlobalEnv[[df]]$logFC.cis*.GlobalEnv[[df]]$logFC.trans>0)]<-"Cis+Trans"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Cis+Trans"]<-4
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<=0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<=0.05)]<-"Compensatory"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Compensatory"]<-5
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Conserved"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Conserved"]<-7
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<=0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)
                            &(.GlobalEnv[[df]]$Category=='')]<-"Ambiguous_total.but.no.cis.or.trans"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.cis.or.trans"]<-6 # large FC
  
  
  .GlobalEnv[[df]]$Category.factor<-as.numeric(.GlobalEnv[[df]]$Category.factor)
  .GlobalEnv[[df]]$Cat.main.factor <- round(.GlobalEnv[[df]]$Category.factor)
  .GlobalEnv[[df]]$Cat.main <- .GlobalEnv[[df]]$Category
  .GlobalEnv[[df]]$Cat.main[.GlobalEnv[[df]]$Cat.main.factor == 1] <- 'Cis related'
  .GlobalEnv[[df]]$Cat.main[.GlobalEnv[[df]]$Cat.main.factor == 2] <- 'Trans related'
  .GlobalEnv[[df]]$Cat.main[.GlobalEnv[[df]]$Cat.main.factor == 2] <- 'Trans related'
  .GlobalEnv[[df]]$Cat.main[.GlobalEnv[[df]]$Cat.main.factor == 6] <- 'Ambiguous'
  .GlobalEnv[[df]]$sub.cat <- .GlobalEnv[[df]]$Category
  .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==1.1] <- "Cis*"
  .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==1.2] <- "Cis**"
  .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==2.1] <- "Trans*"
  .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==2.2] <- "Trans**"
  .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==6] <- "Ambiguous"
  .GlobalEnv[[df]]$dmel_gene<-factor(x = .GlobalEnv[[df]]$dmel_gene,levels = .GlobalEnv[[df]]$dmel_gene[order(.GlobalEnv[[df]]$Category.factor,decreasing = F)])
  
}


write.csv(x = ana.ctrl, file = "ASDE.result.fabo.ctrl.csv",row.names = F)
write.csv(x = ana.wasp, file = "ASDE.result.fabo.wasp.csv",row.names = F)
write.csv(x = ana.reac, file = "ASDE.result.fabo.reac.csv",row.names = F)



#########################################
#####                               #####
#####   DE wihtin sim/sec    #####
#####                               #####
#########################################
### set contrast:
con.DE <- makeContrasts(
  con.DE.sim = SIM.wasp.dsim-SIM.ctrl.dsim,
  con.DE.sec = SEC.wasp.dsec-SEC.ctrl.dsec,
  con.DE.hyb = (HYB.wasp.dsim+HYB.wasp.dsec)-(HYB.ctrl.dsim+HYB.ctrl.dsec),
  levels = design
)
### QL F-test:
for (con in colnames(con.DE)) {
  .GlobalEnv[[paste0('QLFT_DE.',substr(x = con,start = 8,stop = nchar(con)))]] <- glmQLFTest(fit,contrast = con.DE[,con])
  .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]] <- topTags(object = .GlobalEnv[[paste0('QLFT_DE.',substr(x = con,start = 8,stop = nchar(con)))]],n = Inf)$table
  .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]<-.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]%>%
    mutate(dmel_gene = rownames(.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]))
}

### arrnage results:
rownames(df.result_DE.sim)<-c()
rownames(df.result_DE.sec)<-c()
rownames(df.result_DE.hyb)<-c()
ana.DE <- bind_rows(df.result_DE.sim%>%mutate(species='*D. simulans*'),
                    df.result_DE.sec%>%mutate(species='*D. sechellia*'),
                    df.result_DE.hyb%>%mutate(species='Hybrid'))%>%mutate(tissue='Fat Body')



write.csv(x = ana.DE, file = "ASDE.result.fabo.response.within.SIM.SEC.csv",row.names = F)



