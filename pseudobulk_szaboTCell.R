rm(list = ls())

library(Seurat)
library(tidyverse)
library(SCPA)
library(msigdbr)
library(magrittr)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(muscat)

library(scater)
library(scran)



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


saveDir <- file.path(dirPath, "results_pseudobulk")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

list.files(filesPath)

seurat <- read_rds(file.path(filesPath, "szabo_t_cell.rds"))
Idents(seurat)
levels(seurat)

rest <- subset(seurat, subset = stimulation == "none")

# p1 <- DimPlot(rest, reduction = "umap", group.by = "tissue") + NoLegend()
p2 <- RidgePlot(rest, features = "CD3E", group.by = "tissue")
p3 <- RidgePlot(rest, features = "TRAC", group.by = "tissue")
# p1 + 
  p2 + p3

# pseudobulk analysis


sce <- as.SingleCellExperiment(rest)
rm(seurat)
rm(rest)
colData(sce)

p1 <- plotExpression(sce, features = "DEFA1", x = "ident") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                   hjust = 1))
p2 <- plotPCA(sce, colour_by = "ident")
p1 + p2

top <- getTopHVGs(sce, n = 2000)
"DEFA1" %in% top

sce <- prepSCE(sce,
               kid = "fine",
               sid = "tissue",
               gid = "tissue",
               drop = FALSE)
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

t(table(sce$cluster_id, sce$sample_id))
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)
(pb_mds <- pbMDS(pb))
