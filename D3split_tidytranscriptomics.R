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
library(fgsea)
library(hciR)
library(IHW)

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

saveDir <- file.path(dirPath, "results_D3")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

plotDir <- file.path(saveDir, "plots")
dir.create(plotDir)
########################### D3 Hyp ###########################################
########################## prep for data analysis #############################

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

#D3 samples only
grep("D3", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("D3", colnames(rowname_cts))]
rowname_cts["DEFA1", ]

#Hyp samples only
grep("Hyp", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("Hyp", colnames(rowname_cts))]
rowname_cts["DEFA1", ]


# remove DA samples
grep("DA", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]

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

# D3 samples only
grep("D3", rownames(col.data))
col.data <- col.data[grep("D3", rownames(col.data)), ]
col.data

# Hyp samples only
grep("Hyp", rownames(col.data))
col.data <- col.data[grep("Hyp", rownames(col.data)), ]
col.data

# remove DA samples
grep("DA", rownames(col.data))
col.data <- col.data[-grep("DA", rownames(col.data)), ]
col.data

# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))
col.data$Condition <- case_when(col.data$Condition == "Vivo" ~ "Rest",
                      col.data$Condition == "Hyp" ~ "Hypo",
                      col.data$Condition == "Norm" ~ "Norm")
col.data
colnames(col.data)[2] <- "Tissue"
col.data <- col.data %>% mutate(Tissue = str_remove(Tissue, "L"))
col.data

# prep DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Tissue)
design(dds) <- formula(~ Tissue + ID)
dds

########################### tidytranscriptomics ###############################
# coerce DESeqDataSet to RangedSummarizedExperiment
rse <- as(dds, "RangedSummarizedExperiment")

counts <- rse %>% tidybulk()
head(counts)
tail(counts)

# remove "L" from "PBL"
counts_format <- counts %>% 
  mutate(Tissue = str_remove(Tissue, "L"))

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
ggsave(file.path(plotDir, "counts_scaled_Hyp.jpeg"), device = "jpeg")

counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = sample, y = abundance + 1, fill = Tissue)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = median(abundance + 1)), colour="red") +
  facet_wrap(~source) +
  scale_y_log10() +
  theme_bw()
ggsave(file.path(plotDir, "counts_scaled_boxplot_Hyp.jpeg"), device = "jpeg")

counts_scaled %>% group_by(sample) %>%
  summarise(total_reads=sum(counts))

ggplot(counts_scaled, mapping = aes(x = sample, weight = counts, fill = sample)) +
  geom_bar() +
  theme(axis.text.x = element_blank())
ggsave(file.path(plotDir, "count_reads_per_sample_Hyp.jpeg"), device = "jpeg")

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)
counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 100)

counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = Tissue, shape = Condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()
ggsave(file.path(plotDir, "PCA_top100_Hyp.jpeg"), device = "jpeg")
# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- counts_scaled %>% reduce_dimensions(method = "MDS", .dims = 6, top = 100)

tt_mds %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = Tissue, shape = Condition)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.7, type = "norm") +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()
ggsave(file.path(plotDir, "MDS_top100_Hyp.jpeg"), device = "jpeg")


row_labels <- c("DEFA1", "DEFA3", "DEFA4", "ELANE", "CD177", "PRTN3", "MPO", "CXCL12", "FABP4", "PKLR")
name_list <- counts_scaled %>%
  filter(.abundant) %>%
  keep_variable( .abundance = counts_scaled, top = 100) %>%
  as_tibble() %>%
  pull(feature)
name_list <- unique(name_list)
name_list
name_list <- name_list %>% 
  as_tibble() %>%
  mutate(name_list = ifelse(value %in% row_labels, as.character(value), "")) %>%
  pull(name_list)


# hm <- counts_scaled %>%
#   
#   # filter lowly abundant
#   filter(.abundant) %>%
#   
#   # extract most variable genes
#   keep_variable( .abundance = counts_scaled, top = 100) %>%
#   
#   as_tibble() %>%
#   
#   mutate(genes = feature) %>%
# 
#   # create heatmap
#     heatmap(
#     .column = sample,
#     .row = genes,
#     .value = counts_scaled,
#     row_names_gp = gpar(fontsize = 7),
#     transform = log1p,
#     palette_value = c("blue", "white", "red"),
#     show_column_names = FALSE,
#     show_row_names = TRUE,
#     column_km = 2,
#     row_km = 3,
#     row_labels = name_list
#   ) %>%
#   add_tile(c(Tissue))
# hm
hm <- counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 100) %>%
  
  as_tibble() %>%
  
  mutate(genes = feature) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    annotation = c(Tissue),
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 2,
    row_km = 5,
    row_labels = name_list
  )
hm

pdf(file = file.path(plotDir, "heatmap_top100_Hyp_names.pdf"))
hm
dev.off()
# contrasts , "ConditionNorm - ConditionRest", "ConditionHypo - ConditionNorm"
counts_de <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + Tissue + ID,
    .contrasts = c("TissueBM - TissuePB"),
    omit_contrast_in_colnames = TRUE
  )


counts_de %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% arrange(desc(logFC)) %>% pull(feature) %>% head(50)

topgenes <-
  counts_de %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(FDR) %>%
  head(10)

topgenes_symbols <- topgenes %>% pull(feature)

volcano <- counts_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano
pdf(file = file.path(plotDir, "volcanoBM-PBHyp_plot.pdf"))
volcano
dev.off()

# informative MA plot
maplot <- counts_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  
  # Plot
  ggplot(aes(x = logCPM, y = logFC, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
maplot
pdf(file = file.path(plotDir, "MAplot.pdf"))
maplot
dev.off()

counts_de %>%
  filter(.abundant) %>%
  pivot_transcript(.transcript = feature) %>%
  arrange(FDR) %>%
  write_csv(file.path(dirPath, "results", "deBM-PBD3Hyp_results_ordered.csv"))

topgenes_symbols <- c("DEFA1", "DEFA4", "ELANE", "CD177", "CXCL12", "PRTN3", "DEFA3", "FABP4", "MPO")

strip_chart <-
  counts_scaled %>%
  
  # extract counts for top differentially expressed genes
  filter(feature %in% topgenes_symbols) %>%
  
  # make stripchart
  ggplot(aes(x = Tissue, y = counts_scaled + 1, fill = Tissue, label = "")) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ feature) +
  scale_y_log10()+
  theme_bw()

strip_chart

pdf(file = file.path(plotDir, "stripchart_12genes_BM-PBD3Hyp.pdf"))
strip_chart
dev.off()

# deconvolve cellularity
counts_de_cibersort <- deconvolve_cellularity(counts_de, action = "get", cores = 1,
                                              method = "llsr", prefix = "cibersort__")

counts_de_cibersort %>%
  pivot_longer(
    names_to = "Cell_type_inferred", 
    values_to = "proportion", 
    names_prefix ="cibersort__", 
    cols = contains("cibersort__")
  ) %>%
  ggplot(aes( x = `Cell_type_inferred`, y = proportion)) +
  geom_boxplot() +
  facet_wrap(~ sample) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio=1/5)

############################# D3 Norm ##########################################
########################## prep for data analysis #############################

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

#D3 samples only
grep("D3", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("D3", colnames(rowname_cts))]
rowname_cts["DEFA1", ]

#Norm samples only
grep("Norm", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("Norm", colnames(rowname_cts))]
rowname_cts["DEFA1", ]


# remove DA samples
grep("DA", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]

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

# D3 samples only
grep("D3", rownames(col.data))
col.data <- col.data[grep("D3", rownames(col.data)), ]
col.data

# Hyp samples only
grep("Norm", rownames(col.data))
col.data <- col.data[grep("Norm", rownames(col.data)), ]
col.data

# remove DA samples
grep("DA", rownames(col.data))
col.data <- col.data[-grep("DA", rownames(col.data)), ]
col.data

# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))
col.data$Condition <- case_when(col.data$Condition == "Vivo" ~ "Rest",
                                col.data$Condition == "Hyp" ~ "Hypo",
                                col.data$Condition == "Norm" ~ "Norm")
col.data
colnames(col.data)[2] <- "Tissue"
col.data <- col.data %>% mutate(Tissue = str_remove(Tissue, "L"))
col.data

# prep DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Tissue)
design(dds) <- formula(~ Tissue + ID)
dds

########################### tidytranscriptomics ###############################
# coerce DESeqDataSet to RangedSummarizedExperiment
rse <- as(dds, "RangedSummarizedExperiment")

counts <- rse %>% tidybulk()
head(counts)
tail(counts)

# remove "L" from "PBL"
counts_format <- counts %>% 
  mutate(Tissue = str_remove(Tissue, "L"))

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
ggsave(file.path(plotDir, "counts_scaled_Norm.jpeg"), device = "jpeg")

counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = sample, y = abundance + 1, fill = Tissue)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = median(abundance + 1)), colour="red") +
  facet_wrap(~source) +
  scale_y_log10() +
  theme_bw()
ggsave(file.path(plotDir, "counts_scaled_boxplot_Norm.jpeg"), device = "jpeg")

counts_scaled %>% group_by(sample) %>%
  summarise(total_reads=sum(counts))

ggplot(counts_scaled, mapping = aes(x = sample, weight = counts, fill = sample)) +
  geom_bar() +
  theme(axis.text.x = element_blank())
ggsave(file.path(plotDir, "count_reads_per_sample_Norm.jpeg"), device = "jpeg")

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)
counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 100)

counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = Tissue, shape = Condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()
ggsave(file.path(plotDir, "PCA_top500_Norm.jpeg"), device = "jpeg")
# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- counts_scaled %>% reduce_dimensions(method = "MDS", .dims = 6, top = 100)

tt_mds %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = Tissue, shape = Condition)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.7, type = "norm") +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()
ggsave(file.path(plotDir, "MDS_top100_Norm.jpeg"), device = "jpeg")


row_labels <- c("DEFA1", "DEFA3", "DEFA4", "ELANE", "CD177", "PRTN3", "MPO", "CXCL12", "FABP4", "PKLR")
name_list <- counts_scaled %>%
  filter(.abundant) %>%
  keep_variable( .abundance = counts_scaled, top = 100) %>%
  as_tibble() %>%
  pull(feature)
name_list <- unique(name_list)
name_list
name_list <- name_list %>% 
  as_tibble() %>%
  mutate(name_list = ifelse(value %in% row_labels, as.character(value), "")) %>%
  pull(name_list)


# hm <- counts_scaled %>%
#   
#   # filter lowly abundant
#   filter(.abundant) %>%
#   
#   # extract most variable genes
#   keep_variable( .abundance = counts_scaled, top = 100) %>%
#   
#   as_tibble() %>%
#   
#   mutate(genes = feature) %>%
# 
#   # create heatmap
#     heatmap(
#     .column = sample,
#     .row = genes,
#     .value = counts_scaled,
#     row_names_gp = gpar(fontsize = 7),
#     transform = log1p,
#     palette_value = c("blue", "white", "red"),
#     show_column_names = FALSE,
#     show_row_names = TRUE,
#     column_km = 2,
#     row_km = 3,
#     row_labels = name_list
#   ) %>%
#   add_tile(c(Tissue))
# hm
hm <- counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 100) %>%
  
  as_tibble() %>%
  
  mutate(genes = feature) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    annotation = c(Tissue),
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 2,
    row_km = 2
  )
hm

pdf(file = file.path(plotDir, "heatmap_top100_Norm.pdf"))
hm
dev.off()
# contrasts , "ConditionNorm - ConditionRest", "ConditionHypo - ConditionNorm"
counts_de <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + Tissue + ID,
    .contrasts = c("TissueBM - TissuePB"),
    omit_contrast_in_colnames = TRUE
  )


counts_de %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% arrange(desc(logFC)) %>% pull(feature) %>% head(50)

topgenes <-
  counts_de %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(FDR) %>%
  head(10)

topgenes_symbols <- topgenes %>% pull(feature)

volcano <- counts_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano
pdf(file = file.path(plotDir, "volcanoBM-PBNorm_plot.pdf"))
volcano
dev.off()

# informative MA plot
maplot <- counts_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  
  # Plot
  ggplot(aes(x = logCPM, y = logFC, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
maplot
pdf(file = file.path(plotDir, "MAplot.pdf"))
maplot
dev.off()

counts_de %>%
  filter(.abundant) %>%
  pivot_transcript(.transcript = feature) %>%
  arrange(FDR) %>%
  write_csv(file.path(dirPath, "results", "deBM-PBD3Norm_results_ordered.csv"))

topgenes_symbols <- c("DEFA1", "DEFA4", "ELANE", "CD177", "CXCL12", "PRTN3", "DEFA3", "FABP4", "MPO")

strip_chart <-
  counts_scaled %>%
  
  # extract counts for top differentially expressed genes
  filter(feature %in% topgenes_symbols) %>%
  
  # make stripchart
  ggplot(aes(x = Tissue, y = counts_scaled + 1, fill = Tissue, label = "")) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ feature) +
  scale_y_log10()+
  theme_bw()

strip_chart

pdf(file = file.path(plotDir, "stripchart_12genes_BM-PBD3Norm.pdf"))
strip_chart
dev.off()

# deconvolve cellularity
counts_de_cibersort <- deconvolve_cellularity(counts_de, action = "get", cores = 1,
                                              method = "llsr", prefix = "cibersort__")

counts_de_cibersort %>%
  pivot_longer(
    names_to = "Cell_type_inferred", 
    values_to = "proportion", 
    names_prefix ="cibersort__", 
    cols = contains("cibersort__")
  ) %>%
  ggplot(aes( x = `Cell_type_inferred`, y = proportion)) +
  geom_boxplot() +
  facet_wrap(~ sample) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio=1/5)
