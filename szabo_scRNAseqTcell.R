rm(list = ls())

library(Seurat)
library(tidyverse)
library(SCPA)
library(msigdbr)
library(magrittr)
library(ComplexHeatmap)
library(SingleCellExperiment)



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

plotDir <- file.path(saveDir, "plots")
dir.create(plotDir)

list.files(filesPath)

t_cell <- read_rds(file.path(filesPath, "szabo_t_cell.rds"))
DimPlot(t_cell, split.by = "tis_stim", group.by = "fine", ncol = 4)
# pdf(file = file.path(plotDir, "UMAP_all.pdf"))
# DimPlot(t_cell, split.by = "tis_stim", group.by = "fine", ncol = 4)
# dev.off()

##################### Subset t_cell object #####################################
# no stim only
t_cell <- subset(t_cell, subset = stimulation == "none")

t_cell$neat <- case_when(t_cell$tissue == "bl" ~ "Blood",
                         t_cell$tissue == "bm" ~ "BM",
                         t_cell$tissue == "lung" ~ "Lung",
                         t_cell$tissue == "ln" ~ "LN")


########################### VlnPlots ##########################################

# UMAP plot with expression

features <- c("DEFA1", "DEFA3", "ELANE")
FeaturePlot(t_cell, features, split.by = "tissue", ncol = 4)

# pdf(file = file.path(plotDir, "UMAP_featureplot.pdf"))
# FeaturePlot(t_cell, features, split.by = "neat", ncol = 4)
# dev.off()

# violin plots

p1 <- VlnPlot(t_cell, "DEFA1", pt.size = 0.1, group.by = "tissue", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text()) + 
  NoLegend() +
  ylab("Log1p")

p2 <- VlnPlot(t_cell, "DEFA3", pt.size = 0.1, group.by = "tissue", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text()) + 
  NoLegend() +
  ylab("Log1p")

p3 <- VlnPlot(t_cell, "DEFA4", pt.size = 0.1, group.by = "tissue", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text()) + 
  NoLegend() +
  ylab("Log1p")

p4 <- VlnPlot(t_cell, "ELANE", pt.size = 0.1, group.by = "tissue", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text()) + 
  NoLegend() +
  ylab("Log1p")

p5 <- VlnPlot(t_cell, "PRTN3", pt.size = 0.1, group.by = "tissue", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text()) + 
  NoLegend() +
  ylab("Log1p")

p6 <- VlnPlot(t_cell, "MPO", pt.size = 0.1, group.by = "tissue", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text()) + 
  NoLegend() +
  ylab("Log1p")

patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)

pdf(file = file.path(plotDir, "VlnPlots.pdf"))
patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
dev.off()

########################## SCPA ###############################################
# metab <- read_csv(file.path(saveDir, "gene_sets", "combined_metabolic_pathways.csv"))
pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
# consider only HALLMARK and ANTIMICROBIAL pathways 2/2 laptop memory/CPU
pathways <- msigdbr::msigdbr("Homo sapiens") %>%
  filter(grepl("ANTIMICROBIAL", gs_name, ignore.case = TRUE) |
           grepl("PROSTAGLANDIN", gs_name, ignore.case = TRUE) |
           grepl("HALLMARK", x = gs_name, ignore.case = TRUE)) %>%
  format_pathways()
length(pathways)
###################### extract genes from pathways ############################
ex_pathways <- msigdbr::msigdbr("Homo sapiens") %>%
  filter(grepl("ANTIMICROBIAL", gs_name, ignore.case = TRUE)) %>%
  format_pathways()
length(ex_pathways)
is.list(ex_pathways)
# function to extarct genes from pathways
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

geneList <- extract_genes(ex_pathways)

########################## prep Seurat Object for SCPA #########################
cell_types <- unique(t_cell$fine)

split_tissue <- SplitObject(t_cell, split.by = "tissue")
# free some memory
rm(t_cell)
############################# SCPA all cells ###################################
blood <- seurat_extract(split_tissue$bl, 
                        meta1 = "stimulation", value_meta1 = "none")

bm <- seurat_extract(split_tissue$bm, 
                     meta1 = "stimulation", value_meta1 = "none")

ln <- seurat_extract(split_tissue$ln,
                     meta1 = "stimulation", value_meta1 = "none")

lung <- seurat_extract(split_tissue$lung,
                       meta1 = "stimulation", value_meta1 = "none")
# compare all tissues to BM
bm_bl <- compare_pathways(list(bm, blood), pathways)
bm_bl %>% arrange(desc(qval)) %>% as_tibble()
plot_rank(bm_bl, pathway = c("antimicrobial", "apoptosis"))
bm_ln <- compare_pathways(list(bm, ln), pathways)
bm_ln %>% arrange(desc(qval)) %>% as_tibble()
plot_rank(bm_ln, pathway = c("HALLMARK_HYPOXIA", 
                             "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED"))
pdf(file = file.path(plotDir, "SCPA_BMLN.pdf"))
plot_rank(bm_ln, pathway = c("HALLMARK_HYPOXIA", 
                             "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED"))
dev.off()

bm_lung <- compare_pathways(list(bm, lung), pathways)
bm_lung %>% arrange(desc(qval)) %>% as_tibble()
plot_rank(bm_lung, pathway = c("antimicrobial"))

# compare multiple tissues to BM
bm_all <- compare_pathways(list(bm, blood, ln, lung), pathways)
bm_all %>% arrange(desc(qval)) %>% as_tibble()
plot_rank(bm_all, pathway = c("antimicrobial", "hypoxia"))

pdf(file = file.path(plotDir, "SCPA_all.pdf"))
plot_rank(bm_all, pathway = c("antimicrobial", "hypoxia"))
dev.off()


# metabolic pathways only
bm_metab <- compare_pathways(list(bm, blood, ln, lung), 
                             pathway = file.path(saveDir, "gene_sets", "combined_metabolic_pathways.csv"))
bm_metab %>% arrange(desc(qval)) %>% as_tibble()
plot <- plot_rank(bm_metab, pathway = c("HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                                        "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
                                        "REACTOME_PYRUVATE_METABOLISM_AND"))
plot
pdf(file = file.path(plotDir, "SCPA_all_metab.pdf"))
plot
dev.off()

bmln_metab <- compare_pathways(list(bm, blood), 
                             pathway = file.path(saveDir, "gene_sets", "combined_metabolic_pathways.csv"))
bmln_metab %>% arrange(desc(qval)) %>% as_tibble()


plot_rank(bm_bl, pathway = "antimicrobial")
plot_rank(bm_ln, pathway = "antimicrobial")
plot_rank(bm_lung, pathway = "antimicrobial")
colnames(bm_ln)
bm_bl %>% arrange(desc(qval)) %>% select(Pathway, qval, adjPval) %>% as_tibble()
bm_ln %>% arrange(desc(qval)) %>% select(Pathway, qval, adjPval) %>% as_tibble()
bm_lung %>% arrange(desc(qval)) %>% select(Pathway, qval, adjPval) %>% as_tibble()

##################### SCPA by tissue and T cell subset #########################
#create empty lists to store results from the for loop
bl_bm <- list(); bl_ln <- list(); bl_lung <- list()
bm_all <- list()
for (i in cell_types) {
  
  # extract expression data using `seurat_extract` based on tissue, cell_type ("fine"), and stimulation ("none")
  blood <- seurat_extract(split_tissue$bl, 
                          meta1 = "fine", value_meta1 = i,
                          meta2 = "stimulation", value_meta2 = "none")
  
  bm <- seurat_extract(split_tissue$bm, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "stimulation", value_meta2 = "none")
  
  ln <- seurat_extract(split_tissue$ln,
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "stimulation", value_meta2 = "none")

  lung <- seurat_extract(split_tissue$lung,
                         meta1 = "fine", value_meta1 = i,
                         meta2 = "stimulation", value_meta2 = "none")

  # compare all tissues to blood
  print(paste("comparing", i))
  bl_bm[[i]] <- compare_pathways(list(blood, bm), pathways)
  bl_ln[[i]] <- compare_pathways(list(blood, ln), pathways)
  bl_lung[[i]] <- compare_pathways(list(blood, lung), pathways)
  bm_all[[i]] <- compare_pathways(list(bm, blood, ln, lung), pathways)
  
}

get_qvals <- function(scpa_out, name) {
  
  df <- list()
  for (i in names(scpa_out)) {
    df[[i]] <- scpa_out[[i]] %>%
      select(Pathway, qval)
  }
  
  col_names <- names(df)
  for (i in 1:length(df)) {
    df[[i]] <- magrittr::set_colnames(df[[i]], c("Pathway", paste(name, col_names[[i]], sep = "_")))
  }
  
  return(df)
  
}

scpa_results <- Reduce(full_join, c(get_qvals(bl_bm, "bm"),
                                    get_qvals(bl_ln, "ln"),
                                    get_qvals(bl_lung, "lung")))

pway_rest <- scpa_results %>%
  column_to_rownames("Pathway") %>%
  set_colnames(paste("rest", colnames(.), sep = "_")) %>%
  rownames_to_column("Pathway")


all_data <- pway_rest %>%
  column_to_rownames("Pathway")

stim <- all_data %>%
  colnames() %>%
  substr(1, 4) %>%
  str_to_sentence()

tissue <- colnames(all_data) %>%
  substr(6, 9) %>%
  sub(pattern = "_[a-z]", replacement = "") %>%
  gsub(pattern = "bm", replacement = "BM") %>%
  gsub(pattern = "ln", replacement = "LN") %>%
  gsub(pattern = "lung", replacement = "Lung")

lineage <- colnames(all_data) %>%
  str_extract(pattern = "cd[48]") %>%
  str_to_upper() %>%
  replace_na("CD4")

col_an <- HeatmapAnnotation(Stimulation = stim,
                            Tissue = tissue,
                            Lineage = lineage,
                            col = list(Stimulation = c("Rest" = "gray70", "Stim" = "orangered2"),
                                       Tissue = c("BM" = "#cccccc", "LN" = "#be83e6", "Lung" = "#84c476"),
                                       Lineage = c("CD4" = "#4589ff", "CD8" = "#ff6363")),
                            gp = gpar(col = "white", lwd = 0.05),
                            annotation_name_gp = gpar(fontsize = 9),
                            simple_anno_size = unit(3, "mm"))

hm <- all_data %>%
  Heatmap(name = "Qval",
          show_row_names = F, 
          top_annotation = col_an,
          border = T,
          show_row_dend = F,
          show_column_dend = F,
          show_column_names = F)

ht <- draw(hm)

apply(all_data, 1, var) %>%
  data.frame() %>% 
  set_colnames("variation") %>%
  arrange(desc(variation)) %>% 
  rownames_to_column("pathway") %>%
  ggplot(aes(reorder(pathway, variation), variation)) +
  geom_point(shape = 21, cex = 3, fill = "royalblue2", color = 'black', stroke = 0.2) +
  scale_x_discrete(expand = c(0.04, 0.04)) +
  labs(x = "Pathway", y = "Variance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))
plot_rank(scpa_results, pathway = "REACTOME_ANTIMICROBIAL_PEPTIDES")

mydf <- apply(all_data, 1, var) %>%
  data.frame() %>% 
  set_colnames("variation") %>%
  arrange(desc(variation)) %>% 
  rownames_to_column("pathway") %>%
  as_tibble()

rm(split_tissue)






