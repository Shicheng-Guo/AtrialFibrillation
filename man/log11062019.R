#########################################################################################
###   Title : Genome-wide DNA methylation analysis for atrial fibrillation patient
###   Author: Shicheng Guo, Ph.D. Email: Shihcheng.Guo@Gmail.com 
###   Section 1. Function predifinition 
###   Section 2. Data Cleaning
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
###   Section 5. GEO Validation (GSE34639,GSE27895)
###   BIRC10-LC and Windows work well. CHG1 and HPC fail. 
###   11/06/2019
#########################################################################################
BiocManager::install("ChAMP") 
BiocManager::install("doParallel") 
BiocManager::install("benchmarkme") 
benchmarkme::get_ram()
detectCores()

system("cd /home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")
library("ChAMP")
library("doParallel")
Dir="/home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"
set.seed(11)

##############################################################################################################
targets <- read.metharray.sheet(Dir)
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
head(getProbeInfo(manifest))
myNormalRGSet<-preprocessFunnorm(RGSet, nPCs=4, sex = NULL, bgCorr = TRUE,dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,verbose = TRUE)
##############################################################################################################
##############################   Full Automatically Pipeline    ##########################################
##############################################################################################################

testDir="/home/local/MFLDCLIN/guosa/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"
champ.process(directory = testDir)
##############################################################################################################
##############################   Separate Automatically Pipeline    ##########################################
##############################################################################################################
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
seed=110
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=5,method="BMIQ")

# PureG3: 24 case vs 4 control 
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$pureG3,compare.group=c("Case","Control"),arraytype="EPIC")
write.table(myDMP1,file=paste("AtrialFibrillation.",".24Case4Control.pureG3.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$pureG3,method="Bumphunter",compare.group=c("Case","Control"),arraytype="EPIC",minProbes=2,cores=5,maxGap=3000)
write.table(myDMR1,file="AtrialFibrillation.24Case4Control.pureG3.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myGSEA,file="AtrialFibrillation.pureG3.myGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myebayGSEA,file="AtrialFibrillation.pureG3.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")

# Sample_Group: 24 case vs 24 control 
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 0.01)
write.table(myDMP,file=paste("AtrialFibrillation.",".24Case24Control.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=5,maxGap=3000)
write.table(myDMR,file="AtrialFibrillation.24Case24Control.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myGSEA$DMP,file="AtrialFibrillation.DMP.24case24control.myGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myebayGSEA,file="AtrialFibrillation.24case24control.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase <- champ.refbase(beta=myNorm,arraytype="EPIC")
colnames(myRefBase$CorrectedBeta)=myLoad$pd$ID
rownames(myRefBase$CellFraction)=myLoad$pd$ID
write.table(myRefBase$CellFraction,file="myRefBase.CellFraction.txt",sep="\t",quote=F)
myDMPcorrect <- champ.DMP(beta = myRefBase$CorrectedBeta,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 1)
dim(myDMP$Case_to_Control)
P=myDMPcorrect$Case_to_Control[match(rownames(myDMP$Case_to_Control),rownames(myDMPcorrect$Case_to_Control)),]
PP=myDMP$Case_to_Control[match(rownames(P[which(p.adjust(P$P.Value,"bonferroni")<0.05),]),rownames(myDMP$Case_to_Control)),]
write.table(PP,file="Table_1.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 1)
ManhattanmyDMP(myDMP$Case_to_Control)

# Young_Old: 24 case vs 24 control 
myDMP3 <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Young_Old,arraytype="EPIC",adjPVal = 0.1)
write.table(myDMP3,file="AtrialFibrillation.YoungOld.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR3 <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Young_Old,method="Bumphunter",compare.group=c("Old","Young"),arraytype="EPIC",minProbes=2,cores=5,maxGap=3000)
write.table(myDMR3,file="AtrialFibrillation.YoungOld.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMP4 <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Pmr_Enrollment,arraytype="EPIC",adjPVal = 0.1)
write.table(myDMP4,file="AtrialFibrillation.Enrollment.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR4 <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Pmr_Enrollment,method="Bumphunter",compare.group=c("Case","Control"),arraytype="EPIC",minProbes=2,cores=6,maxGap=300)
write.table(myDMR4,file="AtrialFibrillation.Enrollment.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")


myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myebayGSEA,file="AtrialFibrillation.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase1 <- champ.refbase(beta=myNorm,arraytype="EPIC")

myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase1 <- champ.refbase(beta=myNorm,arraytype="EPIC")

myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.CaseControl.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.YoungOld.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.Enrollment.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")



ManhattanmyDMP<-function(myDMP){
  library("qqman")
  library("Haplin")
  SNP=rownames(myDMP)
  CHR=myDMP$CHR
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=myDMP$MAPINFO
  P=myDMP$P.Value
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=0.05/nrow(manhattaninput)
  pdf("manhattan.pdf")
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),lwd=2, suggestiveline=F,genomewideline=FALSE)
  dev.off()
  pdf("qqplot.pdf")
  pQQ(P, nlabs =length(pvalues), conf = 0.95)
  dev.off()
}

 
