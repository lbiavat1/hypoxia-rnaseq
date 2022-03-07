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
library(RColorBrewer)

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

saveDir <- file.path(dirPath, "results_D3-D0")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

plotDir <- file.path(saveDir, "plots")
dir.create(plotDir)
###########################BM D3-D0 ###########################################
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

#BM samples only
grep("BM", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("BM", colnames(rowname_cts))]
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

# BM samples only
grep("BM", rownames(col.data))
col.data <- col.data[grep("BM", rownames(col.data)), ]
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
                              design = ~ ID + Condition)
design(dds) <- formula(~ Condition + ID)
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
# ggsave(file.path(plotDir, "counts_scaled_BM.pdf"), device = "pdf")

counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = sample, y = abundance + 1, fill = Condition)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = median(abundance + 1)), colour="red") +
  facet_wrap(~source) +
  scale_y_log10() +
  theme_bw()
# ggsave(file.path(plotDir, "counts_scaled_boxplot_BM.pdf"), device = "pdf")

counts_scaled %>% group_by(sample) %>%
  summarise(total_reads=sum(counts))

ggplot(counts_scaled, mapping = aes(x = sample, weight = counts, fill = sample)) +
  geom_bar() +
  theme(axis.text.x = element_blank())
# ggsave(file.path(plotDir, "count_reads_per_sample_BM.pdf"), device = "pdf")

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)
counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 100)

counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = Condition, shape = Tissue)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()
# ggsave(file.path(plotDir, "PCA_top500_BM.pdf"), device = "pdf")
# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- counts_scaled %>% reduce_dimensions(method = "MDS", .dims = 6, top = 100)

tt_mds %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = Condition, shape = Tissue)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.7, type = "norm") +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()
# ggsave(file.path(plotDir, "MDS_top100_BM.pdf"), device = "pdf")


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

display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 4, "Paired")

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
    row_names_gp = gpar(fontsize = 4),
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 3,
    column_km_repeats = 100,
    row_km = 4,
    row_km_repeats = 100,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"), font = 1:4)
  ) %>%
  add_tile(c(Condition))
hm

# pdf(file = file.path(plotDir, "heatmap_top100_kmClusters_BM.pdf"))
# hm
# dev.off()

hm <- counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) %>%
  
  as_tibble() %>%
  
  mutate(genes = feature) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 3,
    column_km_repeats = 100,
    row_km = 4,
    row_km_repeats = 100,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"), font = 1:4)
  ) %>%
  add_tile(c(Condition))
hm

# pdf(file = file.path(plotDir, "heatmap_top500_names_BM.pdf"))
# hm
# dev.off()

counts_de_Hyp <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + Condition + ID,
    method = "DESeq2",
    .contrasts = list(c("Condition", "Hypo", "Rest")),
    omit_contrast_in_colnames = TRUE
  )
counts_de_Nor <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + Condition + ID,
    method = "DESeq2",
    .contrasts = list(c("Condition", "Norm", "Rest")),
    omit_contrast_in_colnames = TRUE
  )


hyp <- counts_de_Hyp %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% arrange(desc(stat)) %>% pull(feature) %>% head(50)


nor <- counts_de_Nor %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% arrange(desc(stat)) %>% pull(feature) %>% head(50)

sum(hyp %in% nor)/length(hyp)
################### Activation in Hypoxia ####################################

topgenes_UP <-
  counts_de_Hyp %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(desc(stat)) %>%
  head(10)
topgenes_DN <-
  counts_de_Hyp %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(stat) %>%
  head(10)
topgenes <- full_join(topgenes_DN, topgenes_UP)
topgenes_symbols <- topgenes %>% pull(feature)

volcano <- counts_de_Hyp %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = pvalue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano
# pdf(file = file.path(plotDir, "volcano_BMActivationInHypoxia_plot.pdf"))
# volcano
# dev.off()

# informative MA plot
maplot <- counts_de_Hyp %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  
  # Plot
  ggplot(aes(x = lfcSE, y = log2FoldChange, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
maplot
# pdf(file = file.path(plotDir, "MA_BMActivationInHypoxia_plot.pdf"))
# maplot
# dev.off()

counts_de_Hyp %>%
  filter(.abundant) %>%
  pivot_transcript(.transcript = feature) %>%
  arrange(desc(stat)) %>%
  write_csv(file.path(saveDir, "deBM-D3Hyp-D0_results_ordered.csv"))

################### Activation in Normoxia ####################################

topgenes_UP <-
  counts_de_Nor %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(desc(stat)) %>%
  head(10)
topgenes_DN <-
  counts_de_Nor %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(stat) %>%
  head(10)
topgenes <- full_join(topgenes_DN, topgenes_UP)
topgenes_symbols <- topgenes %>% pull(feature)

volcano <- counts_de_Nor %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = pvalue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano

# pdf(file = file.path(plotDir, "volcano_BMActivationInNormoxia_plot.pdf"))
# volcano
# dev.off()

# informative MA plot
maplot <- counts_de_Nor %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  
  # Plot
  ggplot(aes(x = lfcSE, y = log2FoldChange, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
maplot
# pdf(file = file.path(plotDir, "MA_BMActivationInNormoxia_plot.pdf"))
# maplot
# dev.off()

counts_de_Nor %>%
  filter(.abundant) %>%
  pivot_transcript(.transcript = feature) %>%
  arrange(desc(stat)) %>%
  write_csv(file.path(saveDir, "deBM-D3Nor-D0_results_ordered.csv"))

hh <- counts_de_Hyp %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>%
  arrange(desc(stat)) %>%
  pull(feature)

hh_sig_up <- counts_de_Hyp  %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>% pull(feature)

nn <- counts_de_Nor %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>%
  arrange(desc(stat)) %>%
  pull(feature)

nn_sig_up <- counts_de_Nor %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>% pull(feature)

hh %in% nn_sig_up
hh[!(hh %in% nn_sig_up)] %>% head(10)

nn %in% hh_sig_up
nn[!(nn %in% hh_sig_up)] %>% head(10)

topgenes_symbols <- c(hh[!(hh %in% nn_sig_up)] %>% head(10), nn[!(nn %in% hh_sig_up)] %>% head(10))

strip_chart <-
  counts_scaled %>%
  
  # extract counts for top differentially expressed genes
  filter(feature %in% topgenes_symbols) %>%
  
  # make stripchart
  ggplot(aes(x = Condition, y = counts_scaled + 1, fill = Condition, label = "")) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ feature) +
  scale_y_log10()+
  theme_bw()

strip_chart

# pdf(file = file.path(plotDir, "stripchart_HvsNgenes_BMD3-D0.pdf"))
# strip_chart
# dev.off()

# deconvolve cellularity
counts_de_cibersort <- deconvolve_cellularity(counts_de_Hyp, action = "get", cores = 1,
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

############################# PB D3-D0 ##########################################
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

#PB samples only
grep("PB", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("PB", colnames(rowname_cts))]
rowname_cts["DEFA1", ]

# # remove DA samples
# grep("DA", colnames(rowname_cts))
# rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]
# # remove TP samples
# grep("TP", colnames(rowname_cts))
# rowname_cts <- rowname_cts[, -grep("TP", colnames(rowname_cts))]
# remove MS samples
grep("MS", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("MS", colnames(rowname_cts))]

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

# PB samples only
grep("PB", rownames(col.data))
col.data <- col.data[grep("PB", rownames(col.data)), ]
col.data

# # remove DA samples
# grep("DA", rownames(col.data))
# col.data <- col.data[-grep("DA", rownames(col.data)), ]
# col.data
# # remove TP samples
# grep("TP", rownames(col.data))
# col.data <- col.data[-grep("TP", rownames(col.data)), ]
# col.data

# remove MS samples
grep("MS", rownames(col.data))
col.data <- col.data[-grep("MS", rownames(col.data)), ]
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
                              design = ~ ID + Condition)
design(dds) <- formula(~ Condition + ID)
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
# ggsave(file.path(plotDir, "counts_scaled_PB.pdf"), device = "pdf")

counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = sample, y = abundance + 1, fill = Condition)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = median(abundance + 1)), colour="red") +
  facet_wrap(~source) +
  scale_y_log10() +
  theme_bw()
# ggsave(file.path(plotDir, "counts_scaled_boxplot_PB.pdf"), device = "pdf")

counts_scaled %>% group_by(sample) %>%
  summarise(total_reads=sum(counts))

ggplot(counts_scaled, mapping = aes(x = sample, weight = counts, fill = sample)) +
  geom_bar() +
  theme(axis.text.x = element_blank())
# ggsave(file.path(plotDir, "count_reads_per_sample_PB.pdf"), device = "pdf")

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)
counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 100)

counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = Condition, shape = Tissue)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()
# ggsave(file.path(plotDir, "PCA_top100_PB.pdf"), device = "pdf")
# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- counts_scaled %>% reduce_dimensions(method = "MDS", .dims = 6, top = 500)

tt_mds %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = Condition, shape = Tissue)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.7, type = "norm") +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()
# ggsave(file.path(plotDir, "MDS_top500_PB.pdf"), device = "pdf")


display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 4, "Paired")

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
    row_names_gp = gpar(fontsize = 4),
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 3,
    column_km_repeats = 100,
    row_km = 4,
    row_km_repeats = 100,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"), font = 1:4)
  ) %>%
  add_tile(c(Condition))
hm

# pdf(file = file.path(plotDir, "heatmap_top100_kmClusters_PB.pdf"))
# hm
# dev.off()

hm <- counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) %>%
  
  as_tibble() %>%
  
  mutate(genes = feature) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_km = 3,
    column_km_repeats = 100,
    row_km = 4,
    row_km_repeats = 100,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"), font = 1:4)
  ) %>%
  add_tile(c(Condition))
hm

# pdf(file = file.path(plotDir, "heatmap_top500_PB.pdf"))
# hm
# dev.off()

counts_de_Hyp <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + Condition + ID,
    method = "DESeq2",
    .contrasts = list(c("Condition", "Hypo", "Rest")),
    omit_contrast_in_colnames = TRUE
  )

counts_de_Nor <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + Condition + ID,
    method = "DESeq2",
    .contrasts = list(c("Condition", "Norm", "Rest")),
    omit_contrast_in_colnames = TRUE
  )

hyp <- counts_de_Hyp %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% arrange(desc(stat)) %>% pull(feature) %>% head(50)


nor <- counts_de_Nor %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% arrange(desc(stat)) %>% pull(feature) %>% head(50)

sum(hyp %in% nor)/length(hyp)

################### Activation in Hypoxia ####################################

topgenes_UP <-
  counts_de_Hyp %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(desc(stat)) %>%
  head(10)
topgenes_DN <-
  counts_de_Hyp %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(stat) %>%
  head(10)
topgenes <- full_join(topgenes_DN, topgenes_UP)
topgenes_symbols <- topgenes %>% pull(feature)

volcano <- counts_de_Hyp %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = pvalue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano
# pdf(file = file.path(plotDir, "volcano_PBActivationInHypoxia_plot.pdf"))
# volcano
# dev.off()

# informative MA plot
maplot <- counts_de_Hyp %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  
  # Plot
  ggplot(aes(x = lfcSE, y = log2FoldChange, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
maplot
# pdf(file = file.path(plotDir, "MA_PBMActivationInHypoxia_plot.pdf"))
# maplot
# dev.off()

counts_de_Hyp %>%
  filter(.abundant) %>%
  pivot_transcript(.transcript = feature) %>%
  arrange(desc(stat)) %>%
  write_csv(file.path(saveDir, "dePB-D3Hyp-D0_results_ordered.csv"))

################### Activation in Normoxia ####################################

topgenes_UP <-
  counts_de_Nor %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(desc(stat)) %>%
  head(10)
topgenes_DN <-
  counts_de_Nor %>%
  filter(.abundant) %>%
  pivot_transcript() %>%
  arrange(stat) %>%
  head(10)
topgenes <- full_join(topgenes_DN, topgenes_UP)
topgenes_symbols <- topgenes %>% pull(feature)

volcano <- counts_de_Nor %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = pvalue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano

# pdf(file = file.path(plotDir, "volcano_PBActivationInNormoxia_plot.pdf"))
# volcano
# dev.off()

# informative MA plot
maplot <- counts_de_Nor %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  
  # Plot
  ggplot(aes(x = lfcSE, y = log2FoldChange, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
maplot
# pdf(file = file.path(plotDir, "MA_PBActivationInNormoxia_plot.pdf"))
# maplot
# dev.off()

counts_de_Nor %>%
  filter(.abundant) %>%
  pivot_transcript(.transcript = feature) %>%
  arrange(desc(stat)) %>%
  write_csv(file.path(saveDir, "dePB-D3Nor-D0_results_ordered.csv"))

hh <- counts_de_Hyp %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>%
  arrange(desc(stat)) %>%
  pull(feature)

hh_sig_up <- counts_de_Hyp  %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>% pull(feature)

nn <- counts_de_Nor %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>%
  arrange(desc(stat)) %>%
  pull(feature)

nn_sig_up <- counts_de_Nor %>% filter(.abundant) %>% pivot_transcript(.transcript = feature) %>%
  filter(padj < 0.01) %>%
  filter(abs(log2FoldChange) >=2 ) %>% pull(feature)

hh %in% nn_sig_up
hh[!(hh %in% nn_sig_up)] %>% head(10)

nn %in% hh_sig_up
nn[!(nn %in% hh_sig_up)] %>% head(10)

topgenes_symbols <- c(hh[!(hh %in% nn_sig_up)] %>% head(10), nn[!(nn %in% hh_sig_up)] %>% head(10))

strip_chart <-
  counts_scaled %>%
  
  # extract counts for top differentially expressed genes
  filter(feature %in% topgenes_symbols) %>%
  
  # make stripchart
  ggplot(aes(x = Condition, y = counts_scaled + 1, fill = Condition, label = "")) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ feature) +
  scale_y_log10()+
  theme_bw()

strip_chart

# pdf(file = file.path(plotDir, "stripchart_HvsNgenes_PBD3-D0.pdf"))
# strip_chart
# dev.off()

# deconvolve cellularity
counts_de_cibersort <- deconvolve_cellularity(counts_de_Hyp, action = "get", cores = 1,
                                              method = "cibersort", prefix = "cibersort__")

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
