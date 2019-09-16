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
BiocManager::install("sesame") 
BiocManager::install("doParallel") 
BiocManager::install("benchmarkme") 

benchmarkme::get_ram()
detectCores()

library("sesame") 
library("ChAMP") 
library("sesame") 
library("doParallel") 
library("benchmarkme") 
