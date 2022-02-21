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

list.files(filesPath)

t_cell <- read_rds(file.path(filesPath, "szabo_t_cell.rds"))
DimPlot(t_cell, split.by = "tis_stim", group.by = "fine", ncol = 4)

pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
# pathways <- msigdbr::msigdbr("Homo sapiens") %>%
#   filter(grepl(paste(pws, collapse = "|"), gs_subcat, ignore.case = TRUE) |
#            grepl("HALLMARK", x = gs_name, ignore.case = TRUE)) %>%
#   format_pathways()

# consider only HALLMARK and ANTIMICROBIAL pathways 2/2 laptop memory/CPU
pathways <- msigdbr::msigdbr("Homo sapiens") %>%
  filter(grepl("ANTIMICROBIAL", gs_name, ignore.case = TRUE) |
           grepl("PROSTAGLANDIN", gs_name, ignore.case = TRUE) |
           grepl("HALLMARK", x = gs_name, ignore.case = TRUE)) %>%
  format_pathways()
length(pathways)

cell_types <- unique(t_cell$fine)

split_tissue <- SplitObject(t_cell, split.by = "tissue")
rm(t_cell)

#create empty lists to store results from the for loop
bl_bm <- list(); bl_ln <- list(); bl_lung <- list()
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

mydf <- apply(all_data, 1, var) %>%
  data.frame() %>% 
  set_colnames("variation") %>%
  arrange(desc(variation)) %>% 
  rownames_to_column("pathway") 

rm(split_tissue)
tissue_data <- read_rds(file.path(filesPath, "szabo_t_cell.rds"))
tissue_data$neat <- case_when(tissue_data$tissue == "bl" ~ "Blood",
                              tissue_data$tissue == "bm" ~ "BM",
                              tissue_data$tissue == "lung" ~ "Lung",
                              tissue_data$tissue == "ln" ~ "LN")

plots <- VlnPlot(tissue_data, c("DEFA1", "DEFA3", "CD177", "ELANE", "FABP4"), pt.size = 0, group.by = "neat", ncol = 1, combine = F)

p1 <- VlnPlot(tissue_data, "DEFA1", pt.size = 0, group.by = "neat", ncol = 1) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank()) + 
  NoLegend() +
  ylab("Log1p")

p2 <- VlnPlot(tissue_data, "DEFA3", pt.size = 0, group.by = "neat", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank()) + 
  NoLegend() +
  ylab("Log1p")

p3 <- VlnPlot(tissue_data, "IGHG1", pt.size = 0, group.by = "neat", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank()) + 
  NoLegend() +
  ylab("Log1p")

p4 <- VlnPlot(tissue_data, "CD8A", pt.size = 0, group.by = "neat", ncol = 1) +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank()) + 
  NoLegend() +
  ylab("Log1p")


patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)




