rm(list = ls())

# load libraries
library(DESeq2)
library(hciR)
library(ggplotify)
library(data.table)
library(tidyverse)

# tidyverse-friendly packages
library(tidyHeatmap)
library(tidybulk)
library(ggrepel)
library(plotly)
library(GGally)

# colorblind-friendly packages
library(dittoSeq)

# venn diagram packages
library(VennDiagram)
library(ggVennDiagram)
library(limma)

# gene enrichment analysis packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(fgsea)
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

# where to find gmt files
filePath <- file.path(dirPath, "results")

setwd(savePath)

csvPath <- saveDir
setwd(csvPath)
files <- list.files(path = ".", pattern = ".csv$", full.names = FALSE)
BMvsPB_de <- read_csv(files[1])
BMvsPB_ihw <- read_csv(files[2])

# ============================ Plot Settings ==============================
# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )

############################## Volcano plot #################################
topgenes <-
  BMvsPB_de %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2.0) %>%
  filter(significant) %>%
  arrange(FDR) %>%
  head(20)

topgenes_symbols <- topgenes %>% pull(feature) %>% c(., "FABP4")

plot <- BMvsPB_de %>%
  
  # Subset data
  mutate(FDR = ifelse(is.na(FDR), 1, FDR)) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2.0) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = FDR, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  ggtitle("BM vs PB - Resting") + 
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot

# load msigdbr
library(msigdbr)

hallmark.path <- paste(filePath, "gene_sets", "h.all.v7.1.symbols.gmt", sep = "/")
pathway.HALLMARK <- gmtPathways(hallmark.path)
HALLMARK <- read.gmt(hallmark.path)

hallmark <- msigdbr("Homo sapiens", "H") %>% format_pathways()

c7.path <- paste(filePath, "gene_sets", "c7.all.v7.1.symbols.gmt", sep = "/")
pathway.C7 <- gmtPathways(c7.path)
C7 <- read.gmt(c7.path)


BMvsPB_rnk <- BMvsPB_ihw %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2.0) %>%
  filter(significant) %>%
  dplyr::select(id, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(id) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()

BMvsPB_gH <- GSEA(BMvsPB_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                      TERM2GENE = HALLMARK, by = "fgsea")
BMvsPB_gH@result$ID

BMvsPB_C7 <- GSEA(BMvsPB_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                  TERM2GENE = C7, by = "fgsea")
head(BMvsPB_C7@result$ID)

