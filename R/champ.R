#########################################################################################
###   Title : Genome-wide DNA methylation analysis for atrial fibrillation patient
###   Author: Shicheng Guo, Ph.D. Email: Shihcheng.Guo@Gmail.com 
###   Section 1. Function predifinition 
###   Section 2. Data Cleaning
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
###   Section 5. GEO Validation (GSE34639,GSE27895)
###   BIRC10-LC
#########################################################################################
BiocManager::install("ChAMP") 
BiocManager::install("doParallel") 
BiocManager::install("benchmarkme") 
benchmarkme::get_ram()
detectCores()

library("ChAMP")
library("doParallel")
Dir="/home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"
set.seed(11)

targets <- read.metharray.sheet(Dir)
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
head(getProbeInfo(manifest))
myNormalRGSet<-preprocessFunnorm(RGSet, nPCs=4, sex = NULL, bgCorr = TRUE,dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,verbose = TRUE)

myLoad <- champ.load(Dir,filterBeads=TRUE,arraytype="EPIC")
# EPIC has 411 control probes
pdf("MCaldwell.AMP.EPIC.QC.pdf")
champ.QC()
dev.off()
##########################################################################
pdf("MCaldwell.AMP.EPIC.SVD.pdf")
champ.SVD(beta=myNorm,pd=myLoad$pd)
dev.off()
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))
##########################################################################
# don't use all the cores which will easily be killed by system
detectCores()
seed=sample(seq(1,10000,by=1),1)
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=5)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$pureG3,compare.group=c("Case","Control"),arraytype="EPIC")
write.table(myDMP,file=paste("AtrialFibrillation.",seed,".24Case4Control.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")

myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myDMP,file=paste("AtrialFibrillation.",seed,".CaseControl.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")

myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Young_Old,arraytype="EPIC",adjPVal = 0.1)
write.table(myDMP,file="AtrialFibrillation.YoungOld.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Pmr_Enrollment,arraytype="EPIC",adjPVal = 0.1)
write.table(myDMP,file="AtrialFibrillation.Enrollment.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=6,maxGap=1000)
write.table(myDMR,file="AtrialFibrillation.CaseControl.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Young_Old,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=6,maxGap=1000)
write.table(myDMR,file="AtrialFibrillation.YoungOld.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Pmr_Enrollment,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=6,maxGap=1000)
write.table(myDMR,file="AtrialFibrillation.Enrollment.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.CaseControl.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.YoungOld.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.Enrollment.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myebayGSEA,file="AtrialFibrillation.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase1 <- champ.refbase(beta=myNorm,arraytype="EPIC")
myRefBase2 <- champ.refbase(beta=myNorm,arraytype="450K")

