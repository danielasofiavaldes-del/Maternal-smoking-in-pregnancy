# CellChat v2 ligand-receptor analysis for non-smokers and smokers.
# STB-CTB and CTB-HBC, on spatial proteomics and on snRNA-seq.
# Input: Spatial_Proteomics_Summary_{STB,CTB,HBC}.csv and seu_integrated.rds
# Output: cellchat_*_{ctrl,smk,mapped}.csv

library(CellChat)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)

run_cellchat <- function(object, group.by, min.cells, assay = NULL) {
  if (is.null(assay)) {
    cellchat <- createCellChat(object = object, group.by = group.by)
  } else {
    cellchat <- createCellChat(object, group.by = group.by, assay = assay)
  }
  CellChatDB.use <- subsetDB(
    CellChatDB.human,
    search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
  )
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.1, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  aggregateNet(cellchat)
}

compare_cellchat <- function(cellchat_ctrl, cellchat_smk, label) {
  cellchat_ctrl <- netAnalysis_computeCentrality(cellchat_ctrl)
  cellchat_smk <- netAnalysis_computeCentrality(cellchat_smk)
  object.list <- list(ctrl = cellchat_ctrl, smk = cellchat_smk)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
  fig_flow <- rankNet(cellchat, slot.name = "netP", mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE)
  fig_up <- netVisual_bubble(
    cellchat,
    comparison = c(1, 2),
    max.dataset = 2,
    title.name = "Higher signaling in smoker",
    angle.x = 45,
    remove.isolate = TRUE
  )
  fig_down <- netVisual_bubble(
    cellchat,
    comparison = c(1, 2),
    max.dataset = 1,
    title.name = "Lower signaling in smoker",
    angle.x = 45,
    remove.isolate = TRUE
  )
  cellchat <- identifyOverExpressedGenes(
    cellchat,
    group.dataset = "datasets",
    pos.dataset = "smk",
    features.name = "smk",
    only.pos = FALSE,
    thresh.pc = 0.2,
    thresh.fc = 0.05
  )
  net <- netMappingDEG(cellchat, features.name = "smk")
  communication <- subsetCommunication(cellchat)
  write.csv(communication$ctrl, file = paste0("cellchat_", label, "_ctrl.csv"), row.names = FALSE)
  write.csv(communication$smk, file = paste0("cellchat_", label, "_smk.csv"), row.names = FALSE)
  write.csv(net, file = paste0("cellchat_", label, "_mapped.csv"), row.names = FALSE)
  list(cellchat = cellchat, fig_flow = fig_flow, fig_up = fig_up, fig_down = fig_down)
}

read_proteomics_sce <- function() {
  stb <- read.csv("Spatial_Proteomics_Summary_STB.csv", header = TRUE, check.names = FALSE)
  ctb <- read.csv("Spatial_Proteomics_Summary_CTB.csv", header = TRUE, check.names = FALSE)
  hbc <- read.csv("Spatial_Proteomics_Summary_HBC.csv", header = TRUE, check.names = FALSE)
  protein_data <- merge(stb, ctb, by = "Protein", all = TRUE)
  protein_data <- merge(protein_data, hbc, by = "Protein", all = TRUE)
  protein_data[is.na(protein_data)] <- 0
  metadata_colnames <- c(
    "patient_id", "Smoke_status", "Celltype", "Celltype_sub",
    "raw_file_id", "date", "BATCH", "protein_hits"
  )
  long_names <- colnames(protein_data)[2:length(protein_data)]
  metadata_df <- as.data.frame(do.call(rbind, strsplit(long_names, ", ")))
  colnames(metadata_df) <- metadata_colnames
  counts <- as.matrix(protein_data[2:length(protein_data)])
  colnames(counts) <- sub("^(?:[^,]+,){4}\\s*(\\d+).*", "\\1", long_names)
  rownames(counts) <- protein_data$Protein
  SingleCellExperiment(assays = list(logcounts = counts), colData = metadata_df)
}

if (file.exists("Spatial_Proteomics_Summary_STB.csv")) {
  sce <- read_proteomics_sce()
  sce_tropho <- sce[, colData(sce)$Celltype != "HBC"]
  sce_immune <- sce[, colData(sce)$Celltype != "STB"]
  cellchat_prot_tropho <- compare_cellchat(
    run_cellchat(sce_tropho[, colData(sce_tropho)$Smoke_status == "No_smoke"], "Celltype", min.cells = 2),
    run_cellchat(sce_tropho[, colData(sce_tropho)$Smoke_status == "Smoke"], "Celltype", min.cells = 2),
    "prot_STB_CTB"
  )
  cellchat_prot_immune <- compare_cellchat(
    run_cellchat(sce_immune[, colData(sce_immune)$Smoke_status == "No_smoke"], "Celltype", min.cells = 2),
    run_cellchat(sce_immune[, colData(sce_immune)$Smoke_status == "Smoke"], "Celltype", min.cells = 2),
    "prot_CTB_HBC"
  )
}

if (!exists("seu_int_fil") && file.exists("seu_integrated.rds")) {
  seu_int_fil <- readRDS("seu_integrated.rds")
}
if (exists("seu_int_fil")) {
  Idents(seu_int_fil) <- "cell_type"
  seu_subset <- subset(seu_int_fil, idents = c("CTB", "CTBp", "CTBpf", "HBC", "STB", "STBjuv", "HBCp"))
  Idents(seu_subset) <- "condition"
  seu_cntrl <- subset(seu_subset, idents = "cntrl")
  seu_smk <- subset(seu_subset, idents = "smoker")
  Idents(seu_cntrl) <- "cell_type"
  Idents(seu_smk) <- "cell_type"
  cellchat_rna_tropho <- compare_cellchat(
    run_cellchat(subset(seu_cntrl, idents = c("CTB", "CTBp", "CTBpf", "STB", "STBjuv")), "ident", min.cells = 20, assay = "RNA"),
    run_cellchat(subset(seu_smk, idents = c("CTB", "CTBp", "CTBpf", "STB", "STBjuv")), "ident", min.cells = 20, assay = "RNA"),
    "rna_STB_CTB"
  )
  cellchat_rna_immune <- compare_cellchat(
    run_cellchat(subset(seu_cntrl, idents = c("CTB", "CTBp", "CTBpf", "HBC", "HBCp")), "ident", min.cells = 20, assay = "RNA"),
    run_cellchat(subset(seu_smk, idents = c("CTB", "CTBp", "CTBpf", "HBC", "HBCp")), "ident", min.cells = 20, assay = "RNA"),
    "rna_CTB_HBC"
  )
}
