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

plotDir <- file.path(saveDir, "plots")
dir.create(plotDir)

setwd(savePath)

csvPath <- saveDir
setwd(csvPath)
files <- list.files(path = ".", pattern = ".csv$", full.names = FALSE)
files

BMvsPB_de <- read_csv(files[1])

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

# make it tidy
metab <- read_csv(file.path(saveDir, "gene_sets", "combined_metabolic_pathways.csv"))
colnames(metab) <- c("term", paste0("gene", c(1:739)))
metab <- tidyr::gather(metab, key = genes, value = gene, -term) %>% 
  arrange(term) %>% dplyr::select(term, gene)
metab


# pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
# pathways <- msigdbr::msigdbr("Homo sapiens") %>%
#   filter(grepl("ANTIMICROBIAL", gs_name, ignore.case = TRUE) |
#            grepl("PROSTAGLANDIN", gs_name, ignore.case = TRUE) |
#            grepl("HALLMARK", x = gs_name, ignore.case = TRUE)) %>%
#   dplyr::select(gs_name, gene_symbol) %>%
#   mutate(term = gs_name, gene = gene_symbol) %>%
#   dplyr::select(term, gene)

pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
pathways <- msigdbr::msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pws, collapse = "|"), gs_subcat, ignore.case = TRUE) |
           grepl("HALLMARK", x = gs_name, ignore.case = TRUE)) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(term = gs_name, gene = gene_symbol) %>%
  dplyr::select(term, gene)

# pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
msigdbr_collections() %>% as.data.frame()
pathways_GO <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(term = gs_name, gene = gene_symbol) %>%
  dplyr::select(term, gene)


pathways_GO
length(unique(pathways_GO$term))


BMvsPB_rnk <- BMvsPB_de %>%
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
length(BMvsPB_rnk)


BMvsPB_pathways <- GSEA(BMvsPB_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                      TERM2GENE = pathways, by = "fgsea")
BMvsPB_pathways@result$ID
dotplot(BMvsPB_pathways)
gseaplot2(BMvsPB_pathways, geneSetID = "REACTOME_ANTIMICROBIAL_PEPTIDES")
gseaplot2(BMvsPB_pathways, geneSetID = "REACTOME_NEUTROPHIL_DEGRANULATION")
gseaplot2(BMvsPB_pathways, geneSetID = "REACTOME_INNATE_IMMUNE_SYSTEM")
gseaplot2(BMvsPB_pathways, geneSetID = c("REACTOME_INNATE_IMMUNE_SYSTEM",
                                         "REACTOME_NEUTROPHIL_DEGRANULATION",
                                         "REACTOME_ANTIMICROBIAL_PEPTIDES"))
pdf(file.path(plotDir, "gseaPlot_pathways.pdf"))
gseaplot2(BMvsPB_pathways, geneSetID = c("REACTOME_INNATE_IMMUNE_SYSTEM",
                                         "REACTOME_NEUTROPHIL_DEGRANULATION",
                                         "REACTOME_ANTIMICROBIAL_PEPTIDES"))
dev.off()

BMvsPB_pathways@result$ID
BMvsPB_pathways[grep("ANTIMICROBIAL", BMvsPB_pathways@result$ID)]


BMvsPB_GO <- GSEA(BMvsPB_rnk, exponent = 1, minGSSize = 2, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                  TERM2GENE = pathways_GO, by = "fgsea")

BMvsPB_GO@result$ID
dotplot(BMvsPB_GO)

gseaplot2(BMvsPB_GO, geneSetID = c("GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE" ,
                                   "GOBP_ANTIBACTERIAL_HUMORAL_RESPONSE" ,
                                   "GOBP_KILLING_OF_CELLS_IN_OTHER_ORGANISM_INVOLVED_IN_SYMBIOTIC_INTERACTION" ))

BMvsPB_GO@result["GOBP_DEFENSE_RESPONSE" ,"core_enrichment"]
BMvsPB_GO@result["GOBP_ANTIBACTERIAL_HUMORAL_RESPONSE","core_enrichment"]
unlist(strsplit(BMvsPB_GO@result["GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE","core_enrichment"],
                "/"))

pdf(file.path(plotDir, "gseaPlot_GO.pdf"))
gseaplot2(BMvsPB_GO, geneSetID = c("GOCC_T_CELL_RECEPTOR_COMPLEX",
                                   "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE",
                                   "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE"))
dev.off()


##################### extract genes ###############################
library(SCPA)
# function to extract genes from pathways
extract_genes <- function(pathways){
  genes <- c()
  if(is_tibble(pathways)){
    genes <- pathways$Genes
    genes <- unique(genes)
  }else{
    for(i in 1:length(pathways)){
      genes <- c(genes, pathways[[i]]$Genes)
    }
    genes <- unique(genes)
  }
  return(genes)
}

ex_pathways <- msigdbr::msigdbr("Homo sapiens", category = "C5") %>%
  dplyr::filter(gs_name %in% c("GOCC_T_CELL_RECEPTOR_COMPLEX",
                        "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE",
                        "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE")) %>%
  format_pathways()
gene_list <- extract_genes(ex_pathways)
gene_list


# function to extract relevant genes after GSEA

unlist(strsplit(BMvsPB_GO@result["GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE","core_enrichment"],
                "/"))

