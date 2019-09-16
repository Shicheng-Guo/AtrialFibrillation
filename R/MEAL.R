BiocManager::install("MEAL") 
BiocManager::install("MultiDataSet") 
BiocManager::install("minfiData") 
benchmarkme::get_ram()

library(MEAL)
library(MultiDataSet)
library(minfiData)
library(minfi)
library(ggplot2)

data("MsetEx")
