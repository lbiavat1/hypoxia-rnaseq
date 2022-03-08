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

saveDir <- file.path(dirPath, "results_D3-D0")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

plotDir <- file.path(saveDir, "plots")
dir.create(plotDir)

setwd(savePath)

csvPath <- saveDir
setwd(csvPath)
files <- list.files(path = ".", pattern = ".csv$", full.names = FALSE)
files

BM_HAct <- read_csv(files[1])
BM_NAct <- read_csv(files[2])
PB_HAct <- read_csv(files[3])
PB_NAct <- read_csv(files[4])

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

# ===================== Pathway analysis =======================================

# load msigdbr
library(msigdbr)


######################### prep pathways ########################################
# make it tidy
metab <- read_csv(file.path(saveDir, "gene_sets", "combined_metabolic_pathways.csv"))
colnames(metab) <- c("term", paste0("gene", c(1:739)))
metab <- tidyr::gather(metab, key = genes, value = gene, -term) %>% 
  arrange(term) %>% dplyr::select(term, gene)
metab


hallmark <- msigdbr::msigdbr("Homo sapiens") %>%
  filter(grepl("HALLMARK", x = gs_name, ignore.case = TRUE)) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(term = gs_name, gene = gene_symbol) %>%
  dplyr::select(term, gene)

pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
pathways <- msigdbr::msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pws, collapse = "|"), gs_subcat, ignore.case = TRUE) |
           grepl("HALLMARK", x = gs_name, ignore.case = TRUE)) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(term = gs_name, gene = gene_symbol) %>%
  dplyr::select(term, gene)

# pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
# msigdbr_collections() %>% as.data.frame()
# pathways <- msigdbr(species = "Homo sapiens", category = "C7") %>%
#   dplyr::select(gs_name, gene_symbol) %>%
#   mutate(term = gs_name, gene = gene_symbol) %>%
#   dplyr::select(term, gene)


pathways
length(unique(pathways$term))

########################### prep Ranked Lists ##################################

BM_HAct_rnk <- BM_HAct %>%
  filter(.abundant) %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2.0) %>%
  filter(significant) %>%
  dplyr::select(feature, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(feature) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()

BM_NAct_rnk <- BM_NAct %>%
  filter(.abundant) %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2.0) %>%
  filter(significant) %>%
  dplyr::select(feature, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(feature) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()

PB_HAct_rnk <- PB_HAct %>%
  filter(.abundant) %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2.0) %>%
  filter(significant) %>%
  dplyr::select(feature, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(feature) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()

PB_NAct_rnk <- PB_NAct %>%
  filter(.abundant) %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) >= 2.0) %>%
  filter(significant) %>%
  dplyr::select(feature, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(feature) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()

######################### Pathway analysis ####################################

######################### BM_HAct_rnk ##########################################
# Hallmark

BM_HAct_hallmark <- GSEA(BM_HAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                         TERM2GENE = hallmark, by = "fgsea")
BM_HAct_hallmark@result$ID
dotplot(BM_HAct_hallmark)
pdf(file.path(plotDir, "BM_HAct_dotplot_hallmark.pdf"))
dotplot(BM_HAct_hallmark)
dev.off()

gseaplot2(BM_HAct_hallmark, geneSetID = c("HALLMARK_HYPOXIA",
                                          "HALLMARK_GLYCOLYSIS"  ,
                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
pdf(file.path(plotDir, "BM_HAct_gseaPlot_hallmark.pdf"))
gseaplot2(BM_HAct_hallmark, geneSetID = c("HALLMARK_HYPOXIA",
                                          "HALLMARK_GLYCOLYSIS"  ,
                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
dev.off()

# pathways

BM_HAct_pathways <- GSEA(BM_HAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                      TERM2GENE = pathways, by = "fgsea")
BM_HAct_pathways@result$ID
dotplot(BM_HAct_pathways)

gseaplot2(BM_HAct_pathways, geneSetID = c("WP_TYPE_II_INTERFERON_SIGNALING_IFNG" ,
                                          "HALLMARK_GLYCOLYSIS"  ,
                                          "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"))
pdf(file.path(plotDir, "BM_HAct_gseaPlot_pathways.pdf"))
gseaplot2(BM_HAct_pathways, geneSetID = c("WP_TYPE_II_INTERFERON_SIGNALING_IFNG" ,
                                          "HALLMARK_GLYCOLYSIS"  ,
                                          "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"))
dev.off()

# metab pathways

BM_HAct_metab <- GSEA(BM_HAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                         TERM2GENE = metab, by = "fgsea")
BM_HAct_metab@result$ID

gseaplot2(BM_HAct_metab, geneSetID = c("KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM",
                                          "HALLMARK_GLYCOLYSIS"))

######################### BM_NAct_rnk ##########################################
# Hallmark

BM_NAct_hallmark <- GSEA(BM_NAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                         TERM2GENE = hallmark, by = "fgsea")
BM_NAct_hallmark@result$ID
dotplot(BM_NAct_hallmark)
pdf(file.path(plotDir, "BM_NAct_dotplot_hallmark.pdf"))
dotplot(BM_NAct_hallmark)
dev.off()

gseaplot2(BM_NAct_hallmark, geneSetID = c("HALLMARK_IL2_STAT5_SIGNALING",
                                          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
pdf(file.path(plotDir, "BM_NAct_gseaPlot_hallmark.pdf"))
gseaplot2(BM_NAct_hallmark, geneSetID = c("HALLMARK_IL2_STAT5_SIGNALING",
                                          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
dev.off()

# pathways

BM_NAct_pathways <- GSEA(BM_NAct_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                         TERM2GENE = pathways, by = "fgsea")
BM_NAct_pathways@result$ID
dotplot(BM_NAct_pathways)

gseaplot2(BM_NAct_pathways, geneSetID = c("HALLMARK_INTERFERON_ALPHA_RESPONSE" ,
                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                          "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"))
pdf(file.path(plotDir, "BM_NAct_gseaPlot_pathways.pdf"))
gseaplot2(BM_NAct_pathways, geneSetID = c("HALLMARK_INTERFERON_ALPHA_RESPONSE" ,
                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                          "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"))
dev.off()

# metab pathways

BM_NAct_metab <- GSEA(BM_NAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                      TERM2GENE = metab, by = "fgsea")
# no term enriched?




######################### PB_HAct_rnk ##########################################
# Hallmark

PB_HAct_hallmark <- GSEA(PB_HAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                         TERM2GENE = hallmark, by = "fgsea")
PB_HAct_hallmark@result$ID
dotplot(PB_HAct_hallmark)
pdf(file.path(plotDir, "PB_HAct_dotplot_hallmark.pdf"))
dotplot(PB_HAct_hallmark)
dev.off()
PB_HAct_hallmark <- GSEA(PB_HAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.5,
                         TERM2GENE = hallmark, by = "fgsea")
gseaplot2(PB_HAct_hallmark, geneSetID = c("HALLMARK_HYPOXIA",
                                          "HALLMARK_GLYCOLYSIS"  ,
                                          "HALLMARK_MTORC1_SIGNALING"))
pdf(file.path(plotDir, "PB_HAct_gseaPlot_hypoxia.pdf"))
gseaplot2(PB_HAct_hallmark, geneSetID = c("HALLMARK_HYPOXIA",
                                          "HALLMARK_GLYCOLYSIS"  ,
                                          "HALLMARK_MTORC1_SIGNALING"))
dev.off()

# pathways

PB_HAct_pathways <- GSEA(PB_HAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                         TERM2GENE = pathways, by = "fgsea")
PB_HAct_pathways@result$ID
dotplot(PB_HAct_pathways)

gseaplot2(PB_HAct_pathways, geneSetID = c("HALLMARK_MTORC1_SIGNALING" ,
                                          "REACTOME_DNA_REPLICATION_PRE_INITIATION"  ,
                                          "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"))
pdf(file.path(plotDir, "PB_HAct_gseaPlot_pathways.pdf"))
gseaplot2(PB_HAct_pathways, geneSetID = c("WP_TYPE_II_INTERFERON_SIGNALING_IFNG" ,
                                          "HALLMARK_GLYCOLYSIS"  ,
                                          "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"))
dev.off()

# metab pathways

PB_HAct_metab <- GSEA(PB_HAct_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                      TERM2GENE = metab, by = "fgsea")
PB_HAct_metab@result$ID

gseaplot2(PB_HAct_metab, geneSetID = c("REACTOME_METABOLISM_OF_CARBOHYDRATES",
                                       "HALLMARK_GLYCOLYSIS"))

############################## Venn Diagrams ##################################


