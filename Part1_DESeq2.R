rm(list = ls())

# load libraries 

library(DESeq2)
library(limma)
library(ggplotify)

# tidyverse-friendly
library(tidyverse)
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

# create working directory
workingDir <- "WorkingDirectory"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

saveDir <- file.path(dirPath, "results")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

countFile <- "Counts"

# read count (cts) file
cts <- read_tsv(file.path(filesPath, countFile))
is_tibble(cts)
cts

# prep cts
rowname_cts <- cts
names(rowname_cts)[1] <- "rowname"
rowname_cts <- column_to_rownames(rowname_cts)
colnames(rowname_cts)

# remove D0 samples - analyze D3 only
grep("D0", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("D0", colnames(rowname_cts))]

# load and prep col.data file
coldataFile <- "coldata.txt"

colData <- read_tsv(file.path(filesPath, coldataFile))
is_tibble(colData)
names(colData)[1] <- "rowname"
col.data <- column_to_rownames(colData)

# remove D0 samples
grep("D0", rownames(col.data))
col.data <- col.data[-grep("D0", rownames(col.data)), ]

################### Running DESeq2 ##########################

rowname_cts <- round(rowname_cts)

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Condition)
# setup multifactorial design

# create "group" - ?levels "BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"
dds$group <- factor(paste0(dds$Region, "_", dds$Condition),
                    levels = c("BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"))
design(dds) <- formula(~ID + group)

# Pre-Filtering

dim(dds)
keep <- rowSums( counts(dds) ) >= 25
summary(keep)
dds <- dds[ keep, ]
dim(dds)

# varianceStabilizingTransformation
vsd <- vst(dds, blind = TRUE)
DESeq2::plotPCA(vsd, intgroup = "group", ntop = 500)


# remove outliers (DA) and re-rerun


# remove DA samples
grep("DA", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]
colnames(rowname_cts)


# only integer values
rowname_cts <- round(rowname_cts)


grep("DA", rownames(col.data))
col.data <- col.data[-grep("DA", rownames(col.data)), ]
rownames(col.data)


# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Condition)
# setup multifactorial design

# create "group" - ?levels "BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"
dds$group <- factor(paste0(dds$Region, "_", dds$Condition),
                    levels = c("BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"))
design(dds) <- formula(~ID + group)

# Pre-Filtering

dim(dds)
keep <- rowSums( counts(dds) ) >= 25
summary(keep)
dds <- dds[ keep, ]
dim(dds)

# varianceStabilizingTransformation
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

DESeq2::plotPCA(vsd, intgroup = "group", ntop = 500) + stat_ellipse(type = "norm", level = 0.7)
pcaplot <- DESeq2::plotPCA(vsd, intgroup = "group", ntop = 500) + stat_ellipse(type = "norm", level = 0.7)
class(pcaplot)
# ggsave(file.path(savePath, "PCA_4pts_activated.svg"), plot = pcaplot)
ggsave(file.path(savePath, "PCA_4pts_activated.png"), plot = pcaplot)
