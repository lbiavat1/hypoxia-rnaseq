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
workingDir <- "WorkingDirectory_Szabo"
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

bm_pb <- subset(rest, subset = tissue == c("bm", "bl"))
bm_pb
DimPlot(bm_pb, reduction = "umap", split.by = "tissue")

# bm_pb <- DietSeurat(bm_pb, 
#                     counts = TRUE,
#                     data = TRUE,
#                     scale.data = FALSE,
#                     features = NULL,
#                     assays = NULL,
#                     dimreducs = NULL,
#                     graphs = NULL)
bm_pb

# pseudobulk analysis

seurat <- bm_pb
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts 

metadata <- seurat@meta.data
colnames(metadata)

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- metadata$tissue

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id")]



sce$id <- paste0(sce$stim, sce$ind)
(sce <- prepSCE(sce, 
                kid = "cell", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns

sce <- as.SingleCellExperiment(bm_pb)
rm(seurat)
rm(rest)
rm(bm_pb)
colnames(colData(sce))
# remove undetected genes
dim(sce)
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
sce


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
