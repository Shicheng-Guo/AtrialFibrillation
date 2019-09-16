#########################################################################################
###   Title : Genome-wide DNA methylation analysis for atrial fibrillation patient
###   Author: Shicheng Guo, Ph.D. Email: Shihcheng.Guo@Gmail.com 
###   Section 1. Function predifinition 
###   Section 2. Data Cleaning
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
###   Section 5. GEO Validation (GSE34639,GSE27895)
#########################################################################################
BiocManager::install("ChAMP") 
BiocManager::install("doParallel") 
BiocManager::install("benchmarkme") 
BiocManager::install("FlowSorted.Blood.450k") 

benchmarkme::get_ram()
detectCores()

library("ChAMP")
library("doParallel")
library("FlowSorted.Blood.450k")

Dir="/home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"
set.seed(11)

targets <- read.metharray.sheet(Dir)
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
head(getProbeInfo(manifest))
myNormalRGSet<-preprocessFunnorm(RGSet, nPCs=4, sex = NULL, bgCorr = TRUE,dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,verbose = TRUE)

beta <- getBeta(myNormalRGSet)
phen  <- pData(myNormalRGSet)$Sample_Group

anno<-annotation(myNormalRGSet)

predictedSex <- getSex(myNormalRGSet, cutoff = -2)$predictedSex
anno<-read.table("~/hpc/db/GPL21145-48548_EPIC.txt",head=T,sep="\t")

dmp <- dmpFinder(beta, pheno = phen  , type = "categorical")
dmp.full<-data.frame(dmp,anno[match(rownames(dmp),anno$ID),])
write.table(dmp.full,file="dmp.dmpFinder.minfi.txt",sep="\t",col.names=NA,row.names=T,quote=F)

dmp.sig<-subset(dmp,pval<5.7^-8)
dmp.sig<-data.frame(dmp.sig,anno[match(rownames(dmp.sig),anno$ID),c(1,12,13,7)])
write.table(dmp.sig,file="dmp.dmpFinder.minfi.sig.txt",sep="\t",col.names=NA,row.names=T,quote=F)

designMatrix <- model.matrix(~ phen)
dmr <- bumphunter(myNormalRGSet, design = designMatrix, cutoff = 0.2, B=0, type="Beta")
write.table(dmr$table,file="dmr.bumphunter.minfi.txt",sep="\t",col.names=NA,row.names=T,quote=F)

cellCounts <- estimateCellCounts(RGSet)
write.table(cellCounts,file="cellcounts.estimation.txt",sep="\t",col.names=NA,row.names=T,quote=F)

anno<-read.table("~/hpc/db/GPL21145-48548_EPIC.txt",head=T,sep="\t")

wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.minfi.R -O manhattan.qqplot.minfi.R
Rscript manhattan.qqplot.minfi.R dmp.dmpFinder.minfi.txt
