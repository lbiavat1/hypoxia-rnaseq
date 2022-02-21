rm(list = ls())

if(!require("DESeq2"))
  BiocManager::install("DESeq2")

library(DESeq2)

if(!require("limma"))
  BiocManager::install("limma")

library(limma)

if(!require("ggplotify"))
  BiocManager::install("ggplotify")

library(ggplotify)

tidy_transcriptomics <- TRUE
if(tidy_transcriptomics)
  BiocManager::install(c("tidybulk", "tidyHeatmap", "ggrepel", "plotly", "GGally"))

library(tidyHeatmap)
library(tidybulk)
library(ggrepel)
library(plotly)
library(GGally)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory