# Computational panels for Figures 1-5 and Supplementary Figures 2 and 4.
# Reads seu_integrated.rds, sce_umap.rds, and the result tables from scripts 03-07.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(Augur)
library(slingshot)

if (!exists("seu_int_fil") && file.exists("seu_integrated.rds")) {
  seu_int_fil <- readRDS("seu_integrated.rds")
}
if (!exists("deg_results") && file.exists("deg_pseudobulk.csv")) {
  deg_results <- read.csv("deg_pseudobulk.csv")
}
if (!exists("sce_umap") && file.exists("sce_umap.rds")) {
  sce_umap <- readRDS("sce_umap.rds")
}
if (!exists("da_STB") && file.exists("da_STB.csv")) {
  da_STB <- read.csv("da_STB.csv", check.names = FALSE)
  da_CTB <- read.csv("da_CTB.csv", check.names = FALSE)
  da_HBC <- read.csv("da_HBC.csv", check.names = FALSE)
}
if (!exists("da_stb_d6") && file.exists("da_stb_d6.csv")) {
  da_stb_d6 <- read.csv("da_stb_d6.csv", check.names = FALSE)
}

cluster_colors <- c(
  "CCT" = "#D2042D", "EB" = "#FAA0A0", "CTB" = "#66C2A5", "STB" = "#4682B4",
  "FB" = "#954535", "CTBp" = "#FFAA33", "CTBpf" = "#008080", "MC" = "#F781BF",
  "STBjuv" = "#40E0D0", "PAMM" = "#FF7F00", "HBCp" = "#32CD32", "HBC" = "#E97451",
  "VEC" = "#6A3D9A"
)

paper <- c(
  "PEG10", "TP63", "LGR5", "ERVFRD-1", "ERVV-1", "MKI67", "TOP2A", "CYP11A1", "CYP19A1",
  "CGA", "KISS1", "CSH2", "HLA-G", "NOTUM", "CD14", "F13A1", "LYVE1", "CD163", "HLA-DRA",
  "HLA-B", "HBA1", "ANK1", "COL1A1", "COL6A2", "SOX5", "DLK1", "AGTR1", "GUCY1A2",
  "CD34", "DACH1", "MEOX2", "PECAM1"
)

volcano_celltype <- function(deg, celltype) {
  small <- deg %>% filter(cluster == celltype)
  lab <- small %>% filter(FDR < 0.1) %>% slice_max(abs(logFC), n = 20) %>% pull(gene)
  EnhancedVolcano(
    small,
    lab = small$gene,
    x = "logFC",
    y = "PValue",
    title = paste(celltype, "smoking DEG"),
    selectLab = lab,
    pCutoff = 0.1,
    FCcutoff = 0.25,
    pointSize = 1.0,
    labSize = 3.0,
    labCol = "black",
    drawConnectors = TRUE,
    boxedLabels = FALSE,
    col = c("black", "chartreuse4", "darkslategray", "red3"),
    colAlpha = 0.7
  )
}

## Figure 1B. UMAP of annotated nuclei
if (exists("seu_int_fil")) {
  fig1b <- DimPlot(
    seu_int_fil,
    reduction = "umap",
    group.by = "cell_type",
    pt.size = 1.5,
    raster = TRUE,
    label = FALSE
  ) + scale_colour_manual(values = cluster_colors)

  ## Figure 1C. Canonical markers
  DefaultAssay(seu_int_fil) <- "RNA"
  seu_int_fil <- NormalizeData(seu_int_fil, verbose = FALSE)
  fig1c <- DotPlot(seu_int_fil, features = paper, cluster.idents = FALSE, dot.scale = 5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    coord_flip()

  ## Figure 1E. Augur rank on the integrated object
  expr <- GetAssayData(seu_int_fil, slot = "data", assay = "RNA")
  augur_result <- calculate_auc(
    expr,
    seu_int_fil@meta.data,
    cell_type_col = "cell_type",
    label_col = "condition"
  )
  fig1e <- plot_umap(
    augur_result,
    seu_int_fil,
    mode = "rank",
    reduction = "umap",
    cell_type_col = "cell_type",
    palette = "cividis"
  )

  ## Figure 1F. Nuclei proportions by smoking status and gestational age
  cell_abund_df <- prop.table(table(seu_int_fil$cell_type, seu_int_fil$sample_id), margin = 2) %>%
    as.data.frame() %>%
    rename(cell_type = Var1, sample_id = Var2, abundance = Freq)
  sample_meta <- seu_int_fil@meta.data %>%
    select(sample_id, GA, condition) %>%
    distinct()
  cell_abund_df <- cell_abund_df %>%
    left_join(sample_meta, by = "sample_id") %>%
    mutate(
      GA = as.numeric(as.character(GA)),
      GA_group = case_when(
        GA < "56" ~ "GA1 (4–7w)",
        GA > "56" ~ "GA2 (8–11w)",
        TRUE ~ NA_character_
      )
    )
  pd <- position_dodge(width = 0.7)
  fig1f <- ggplot(cell_abund_df, aes(x = cell_type, y = abundance * 100, fill = condition)) +
    geom_boxplot(position = pd, width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black") +
    geom_point(aes(color = condition), position = pd, size = 1.5) +
    facet_grid(GA_group ~ ., scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c("cntrl" = "grey50", "smoker" = "darkred")) +
    scale_color_manual(values = c("cntrl" = "grey50", "smoker" = "darkred")) +
    labs(x = NULL, y = "Abundance (%)", fill = "Condition", color = "Condition") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

## Figure 1D. DEG counts at 10% FDR, and mean nuclei proportion
if (exists("deg_results")) {
  fig1d_deg <- deg_results %>%
    filter(FDR < 0.1) %>%
    group_by(cluster) %>%
    summarise(up = sum(logFC > 0), down = sum(logFC < 0), .groups = "drop") %>%
    pivot_longer(c(up, down), names_to = "direction", values_to = "n") %>%
    ggplot(aes(x = reorder(cluster, n), y = n, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(up = "#7B2D8E", down = "#3B6EA5")) +
    labs(x = NULL, y = "DEGs at FDR < 0.1")

  ## Figure 2A and 2D. STB and HBC volcano plots
  fig2a <- volcano_celltype(deg_results, "STB")
  fig2d <- volcano_celltype(deg_results, "HBC")

  ## Supplementary Figure 2D. logFC of genes differential in two or more groups
  htmap <- deg_results %>% filter(FDR <= 0.1)
  selected_data <- htmap %>% group_by(gene) %>% filter(n() >= 2) %>% ungroup()
  mat <- selected_data %>% pivot_wider(id_cols = cluster, names_from = gene, values_from = logFC)
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$cluster
  mat <- mat[, -1, drop = FALSE]
  fig_s2d <- pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("navy", "white", "red4"))(20),
    na_col = "gray96",
    main = "DEGs in two or more cell types (FDR <= 0.1)"
  )
}

if (exists("seu_int_fil")) {
  fig1d_nuclei <- cell_abund_df %>%
    group_by(cell_type) %>%
    summarise(mean_pct = mean(abundance) * 100, sem = sd(abundance) / sqrt(n()) * 100, .groups = "drop") %>%
    ggplot(aes(x = reorder(cell_type, mean_pct), y = mean_pct)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_pct - sem, ymax = mean_pct + sem), width = 0.2) +
    coord_flip() +
    labs(x = NULL, y = "Mean nuclei proportion (%)")
}

## Figure 3F to 3H. Trophoblast UMAP, pseudotime, and smoking comparison
if (exists("sce_umap")) {
  umap_df <- as.data.frame(reducedDims(sce_umap)$UMAP)
  umap_df$cell_type <- colData(sce_umap)$cell_type
  fig3f <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = cell_type)) +
    geom_point(size = 0.3) +
    scale_colour_manual(values = cluster_colors) +
    theme_classic()
  pt <- slingPseudotime(sce_umap)
  pt_df <- as.data.frame(pt)
  pt_df$cell <- rownames(pt_df)
  pt_df$condition <- colData(sce_umap)$condition
  pt_df$cell_type <- colData(sce_umap)$cell_type
  fig3g <- pt_df %>%
    pivot_longer(starts_with("Lineage"), names_to = "lineage", values_to = "pseudotime") %>%
    filter(!is.na(pseudotime)) %>%
    ggplot(aes(x = pseudotime, y = lineage, color = cell_type)) +
    geom_point(size = 0.4) +
    scale_colour_manual(values = cluster_colors)
  fig3h <- pt_df %>%
    pivot_longer(starts_with("Lineage"), names_to = "lineage", values_to = "pseudotime") %>%
    filter(!is.na(pseudotime)) %>%
    ggplot(aes(x = condition, y = pseudotime, fill = condition)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ lineage) +
    scale_fill_manual(values = c("cntrl" = "grey50", "smoker" = "darkred"))
}

## Figure 4G and Figure 5D. Counts and STB volcano of spatial proteomics
protein_volcano <- function(da, title) {
  da_plot <- da
  da_plot$log2FoldChange <- da_plot$logFC
  da_plot$baseMean <- da_plot$AveExpr
  da_plot$padj <- da_plot$adj.P.Val
  rownames(da_plot) <- da_plot$Protein
  ggmaplot(
    da_plot,
    main = title,
    fdr = 0.05,
    fc = 0.25,
    genenames = da_plot$Protein,
    size = 2.5,
    alpha = 0.5,
    font.label = c(10, "plain", "black"),
    palette = c("#B31B21", "#1465AC", "darkgray"),
    top = 20,
    select.top.method = c("padj", "fc"),
    xlab = "Log2 mean expression",
    ylab = "Log2 fold change",
    ggtheme = ggplot2::theme_classic()
  )
}

if (exists("da_STB") && exists("da_CTB") && exists("da_HBC")) {
  da_counts <- bind_rows(
    da_STB %>% mutate(celltype = "STB"),
    da_CTB %>% mutate(celltype = "CTB"),
    da_HBC %>% mutate(celltype = "HBC")
  ) %>%
    filter(adj.P.Val < 0.05) %>%
    mutate(direction = ifelse(logFC > 0, "up", "down")) %>%
    count(celltype, direction)
  fig4g <- ggplot(da_counts, aes(x = celltype, y = n, fill = direction)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c(up = "#C51B7D", down = "#2166AC")) +
    labs(x = NULL, y = "Proteins at FDR < 0.05")
  fig5d <- protein_volcano(da_STB, "Dysregulated proteins in microdissected STB")
}

## In-vitro STB volcano, from the bulk proteomics EnhancedVolcano call
if (exists("da_stb_d6")) {
  results2 <- da_stb_d6
  keyvals <- ifelse(results2$adj.P.Val < 0.05, "red4", "grey40")
  keyvals[is.na(keyvals)] <- "grey40"
  fig_invitro <- EnhancedVolcano(
    results2,
    lab = results2$Protein,
    x = "logFC",
    y = "P.Value",
    title = "In-vitro STB proteomic differences after CSE",
    selectLab = results2$Protein[results2$adj.P.Val < 0.05],
    pointSize = 2.5,
    labSize = 3.5,
    colCustom = keyvals
  )
}
