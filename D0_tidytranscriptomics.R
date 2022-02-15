rm(list = ls())


# tidyverse-friendly
library(tidyverse)
library(tidybulk)
library(tidyHeatmap)
library(ggrepel)
library(plotly)
library(GGally)

# other useful libraries
library(DESeq2)
library(limma)
library(edgeR)

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
head(rowname_cts)

#D0 samples only
grep("D0", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("D0", colnames(rowname_cts))]
rowname_cts["DEFA1", ]
rowname_cts["LEP", ]

# only integer values
rowname_cts <- round(rowname_cts)
rowname_cts["DEFA1", ]
rowname_cts["DEFA4", ]
rowname_cts["LEP", ]

# load and prep col.data file
coldataFile <- "coldata.txt"

colData <- read_tsv(file.path(filesPath, coldataFile))
is_tibble(colData)
names(colData)[1] <- "rowname"
col.data <- column_to_rownames(colData)

# remove D0 samples
grep("D0", rownames(col.data))
col.data <- col.data[grep("D0", rownames(col.data)), ]

# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))
col.data$Condition <- "Resting"
colnames(col.data)[2] <- "Tissue"
col.data

# prep DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Tissue)
design(dds) <- formula(~ Tissue + ID)
dds
# coerce DESeqDataSet to RangedSummarizedExperiment
rse <- as(dds, "RangedSummarizedExperiment")

counts <- rse %>% tidybulk()
head(counts)
tail(counts)

# remove "L" from "PBL"
counts_format <- counts %>% 
  mutate(Region = str_remove(Tissue, "L"))

counts_scaled <- counts_format %>%
  identify_abundant(factor_of_interest = Tissue, minimum_counts = 25, minimum_proportion = 0.25) %>%
  scale_abundance(method = "TMM")

counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10() +
  theme_bw()

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)

counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = Tissue, shape = Tissue)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()

# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- counts_scaled %>% reduce_dimensions(method = "MDS", .dims = 6, top = 500)

tt_mds %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = Tissue, shape = Tissue)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.7, type = "norm") +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()

counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract 500 most variable genes
  keep_variable( .abundance = counts_scaled, top = 30) %>%
  
  as_tibble() %>%

  # create heatmap
    heatmap(
    .column = sample,
    .row = feature,
    .value = counts_scaled,
    annotation = c(Tissue),
    transform = log1p
  )
