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
