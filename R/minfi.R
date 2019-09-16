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

predictedSex <- getSex(myNormalRGSet, cutoff = -2)$predictedSex

dmp <- dmpFinder(beta, pheno = phen  , type = "categorical")
designMatrix <- model.matrix(~ phen)
dmr <- bumphunter(myNormalRGSet, design = designMatrix, cutoff = 0.2, B=0, type="Beta")

cellCounts <- estimateCellCounts(RGSet)
