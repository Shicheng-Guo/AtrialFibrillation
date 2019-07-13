#####################################################################
###   Title : Genome-wide DNA methylation analysis for RA patient
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Section 1. Function predifinition 
###   Section 2. Data Cleaning
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
###   Section 5. GEO Validation (GSE34639,GSE27895)
#####################################################################
####################################################################################################################################
### Section 1. function predifinition 
####################################################################################################################################
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("minfiData")
BiocManager::install("missMethyl")
BiocManager::install("minfiData")
BiocManager::install("Gviz")
BiocManager::install("DMRcate")
library("ggplot2")
require("minfi")
library("knitr")
library("limma")
library("minfi")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kmanifest")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("IlluminaHumanMethylationEPICmanifest")
library("RColorBrewer")
library("missMethyl")
library("minfiData")
library("Gviz")
library("DMRcate")
library("stringr")
library(minfi)
library(minfiData)
library(sva)
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")


RINfun=function(yorig){
  yranks=rank(yorig)
  tempp=(yranks-.5)/(length(yranks))
  return(qnorm(tempp))
}
RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
}   
############################################################################################
### Section 2. read the idat
############################################################################################
baseDir="C:\\Users\\shg047\\Documents\\GitHub\\AtrialFibrillation\\idat"
baseDir="/mnt/bigdata/Genetic/Projects/shg047/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"

setwd(baseDir)
list.files()
dataDirectory <- baseDir
list.files(dataDirectory, recursive = TRUE)
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(base = baseDir, targets = targets)
pdf("Figure_S1.pdf")
densityPlot(RGset,xlim=c(0,1),sampGroups = RGset$Sample_Group,main = "Beta", xlab = "Beta",cex=0.1)
detP <- detectionP(RGset)
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP), col=pal[factor(targets$Case_Control)], las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Case_Control)), fill=pal,bg="white")
dev.off()

qcReport(RGset, sampNames=targets$ID, sampGroups=targets$Case_Control,pdf="Figure_S2.qcReport.pdf")
head(targets)
keep <- colMeans(detP) < 0.05
rgSet <- RGset[,keep]
rgSet
targets <- targets[keep,]
targets[,1:5]
mSetSq <- preprocessQuantile(rgSet) 
mSetRaw <- preprocessRaw(rgSet)

pdf("Figure_S4.pdf")
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Case_Control,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Case_Control)), text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Case_Control,main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Case_Control)), text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("Figure_S5.pdf")
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)],pch=16,cex=1.5)
legend("top", legend=levels(factor(targets$Case_Control)), text.col=pal,bg="white", cex=0.7,pch=16,col=pal)
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Young_Old)],pch=16,cex=1.5)
legend("top", legend=levels(factor(targets$Young_Old)), text.col=pal,bg="white", cex=0.7,pch=16,col=pal)
dev.off()

pdf("Figure_S6.pdf")
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)], dim=c(1,3),pch=16,cex=1.5)
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,cex=0.7, bg="white",pch=16,col=pal)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Case_Control)], dim=c(2,3),pch=16,cex=1.5)
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,cex=0.7, bg="white",pch=16,col=pal)

plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)], dim=c(3,4),pch=16,cex=1.5)
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,cex=0.7, bg="white",pch=16,col=pal)
dev.off()

mSetSqFlt<-mSetSq
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

head(bVals[,1:5])
head(mVals[,1:5])

pdf("Figure_S6.pdf")
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Case_Control, main="Beta values", legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Case_Control)),text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Case_Control, main="M-values",legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Case_Control)), text.col=brewer.pal(8,"Dark2"))
dev.off()

ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
####################################################################################################################################
### Section 3. DMS and DMR analysis
####################################################################################################################################
MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE,normalize = "controls", reference = 2)
mdsPlot(MSet.norm, numPositions = 1000, sampGroups = MSet.norm$Sample_Group, sampNames = MSet.norm$Sample_Name)
mset <- MSet.norm
M <- getM(mset, type = "beta", betaThreshold = 0.001)
dmp <- dmpFinder(M, pheno=mset$Case_Control, type="categorical",qCutoff=0.05)

write.table(dmp,file="AtrialFibrillation.DMR.txt",sep="\t",col.names = NA,row.names = T,quote=F)

pdf("Figure_S6_dmp.pdf")
plotCpg(mset, cpg=rownames(dmp)[1], pheno=mset$Case_Control)
dev.off()
mset <- mset[rownames(dmp),]
mse <- mapToGenome(mset)
rowData(mse)
mcols(rowData(mse)) <- cbind(mcols(rowData(mse)), dmp)

mSetSq <- preprocessQuantile(rgSet) 
mSetRaw <- preprocessRaw(rgSet)

pdf("Figure_S7.pdf")
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)])
legend("top", legend=levels(factor(targets$Case_Control)), text.col=pal,bg="white", cex=0.7)
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,bg="white", cex=0.7)
dev.off()

pdf("Figure_S8.pdf")
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Case_Control)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Case_Control)), text.col=pal, cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Case_Control)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Case_Control)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Case_Control)), text.col=pal, cex=0.7, bg="white")
dev.off()
