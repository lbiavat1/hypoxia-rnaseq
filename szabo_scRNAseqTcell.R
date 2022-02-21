rm(list = ls())

library(Seurat)
library(tidyverse)
library(SCPA)



# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "WorkingDirectory_D0"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

saveDir <- file.path(dirPath, "results")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

list.files(filesPath)

t_cell <- read_rds(file.path(filesPath, "szabo_t_cell.rds"))
DimPlot(t_cell, split.by = "tis_stim", group.by = "fine", ncol = 4)
