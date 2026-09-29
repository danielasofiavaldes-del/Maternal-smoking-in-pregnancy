# Nucleus QC, doublet scoring, and SCTransform integration on the CellBender counts.
# Input: *_cellbender.rds from 00b_cellbender_remove_background.sh, and cell_cycle.rda
# Output: seu_integrated.rds

library(Seurat)
library(plyr)
library(tidyverse)
library(glmGamPoi)
library(sctransform)
library(DoubletFinder)
library(fields)
library(KernSmooth)
library(ROCR)

options(future.globals.maxSize = +Inf)

# One Seurat object per sample. Barcodes are prefixed with sample_id so they stay unique after merging.
read_single_cell_data <- function(input_data_pattern = "_cellbender.rds") {
  list <- list.files(pattern = "_cellbender\\.rds$", full.names = TRUE)
  list_names <- gsub("_cellbender.rds", "", basename(list))
  data_list <- lapply(list, readRDS)
  data_list <- lapply(seq_along(list_names), function(i) {
    seu <- data_list[[i]]
    seu$sample_id <- list_names[i]
    seu <- RenameCells(seu, add.cell.id = list_names[i])
    seu
  })
  data_list
}

Preprocessing <- function(list) {
  lapply(list, function(seu) {
    seu[["log10GenesPerUMI"]] <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
    seu[["percent_mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu[["percent_ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
    seu[["percent_hb"]] <- PercentageFeatureSet(seu, pattern = "^HB[^(P)]")
    seu[["mitoRatio"]] <- seu@meta.data$percent_mt / 100
    seu[["Percent.Largest.Gene"]] <- apply(
      seu@assays$RNA@counts, 2, function(x) (100 * max(x)) / sum(x)
    )
    seu
  })
}

# percent_mt is on the 0-100 scale: < 0.5 means less than 0.5% mitochondrial counts.
# Upper gene count is 8500. Largest gene under 30.
filter <- function(list) {
  lapply(list, function(seu) {
    seu <- subset(seu, subset = (nFeature_RNA >= 200 & nFeature_RNA <= 8500 & log10GenesPerUMI > 0.80 & percent_mt < 0.5 & Percent.Largest.Gene < 30))
    seu <- seu[!grepl("MALAT1", rownames(seu)), ]
    counts <- GetAssayData(seu, slot = "counts")
    nonzero <- counts > 0
    keep_genes <- Matrix::rowSums(nonzero) >= 10
    filtered_counts <- counts[keep_genes, ]
    CreateSeuratObject(filtered_counts, meta.data = seu@meta.data)
  })
}

# DoubletFinder per sample, for QC.
# Clustering for the homotypic doublet estimate is resolution 0.1.
# nExp is the pK estimate adjusted by the homotypic proportion.
# The calls are written to a table and kept on the object. Nuclei are not removed.
# Sequencing depth is not regressed.
FindDoublets <- function(list) {
  lapply(list, function(seu) {
    seu <- SCTransform(seu, method = "glmGamPoi")
    seu <- RunPCA(seu, dims = 1:50)
    stdv <- seu[["pca"]]@stdev
    percent.stdv <- (stdv / sum(stdv)) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
    co2 <- sort(which((percent.stdv[1:(length(percent.stdv) - 1)] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = TRUE)[1] + 1
    min.pc <- min(co1, co2)
    seu <- RunUMAP(seu, dims = 1:min.pc, reduction = "pca")
    seu <- FindNeighbors(seu, dims = 1:min.pc)
    seu <- FindClusters(seu, resolution = 0.1)
    sweep.res.list <- paramSweep_v3(seu, PCs = 1:min.pc, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
    optimal.pk <- bcmvn.max$pK
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    nExp_poi <- round(optimal.pk * nrow(seu@meta.data))
    nExp_poi_adj <- round(nExp_poi * (1 - modelHomotypic(seu@meta.data$seurat_clusters)))
    seu <- doubletFinder_v3(seu, PCs = 1:min.pc, pK = optimal.pk, nExp = nExp_poi_adj, sct = TRUE)
    DF.name <- colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]
    write.csv(seu@meta.data[, DF.name, drop = FALSE], file = paste0("doublets_", seu$sample_id[1], ".csv"))
    seu
  })
}

# Integrate on 4000 variable genes.
# Cell-cycle scores are stored. Mitochondrial counts and sequencing depth are not regressed.
# The integrated assay is for clustering. Later differential expression uses RNA counts.
integrate_SCT_single_cell_data <- function(seul, cycle_genes) {
  load(cycle_genes)
  for (i in seq_along(seul)) {
    seul[[i]] <- CellCycleScoring(seul[[i]], g2m.features = g2m_genes, s.features = s_genes)
    seul[[i]] <- SCTransform(seul[[i]], method = "glmGamPoi")
    DefaultAssay(seul[[i]]) <- "SCT"
  }
  features <- SelectIntegrationFeatures(object.list = seul, nfeatures = 4000)
  seul <- PrepSCTIntegration(object.list = seul, anchor.features = features)
  integ_anchors <- FindIntegrationAnchors(
    object.list = seul,
    anchor.features = features,
    normalization.method = "SCT"
  )
  seu_int <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
  DefaultAssay(seu_int) <- "integrated"
  seu_int
}

sample_ids <- c("009", "010", "043", "036", "035", "042", "040", "011", "012", "039", "041", "037", "034", "038")

data_list <- read_single_cell_data()
data_list <- data_list[order(vapply(data_list, function(seu) seu$sample_id[1], character(1)))]
data_list <- Preprocessing(data_list)
data_list <- filter(data_list)
# Samples 038, 040 and 043 did not pass QC and are dropped before doublet detection.
data_list <- data_list[!vapply(data_list, function(seu) seu$sample_id[1] %in% c("038", "040", "043"), logical(1))]
datal_filtered_nd <- FindDoublets(data_list)

seu_int <- integrate_SCT_single_cell_data(datal_filtered_nd, cycle_genes = "cell_cycle.rda")

seu_int$GA <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("68", "70", "75", "49", "52", "61", "65", "74", "76", "35", "43", "58", "60", "64"))
seu_int$condition <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("smoker", "smoker", "smoker", "smoker", "smoker", "smoker", "smoker", "cntrl", "cntrl", "cntrl", "cntrl", "cntrl", "cntrl", "cntrl"))
seu_int$gender <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("female", "female", "male", "female", "male", "male", "female", "female", "female", "female", "female", "male", "male", "female"))
seu_int$batch <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("1", "1", "4", "3", "2", "4", "3", "1", "1", "3", "3", "3", "2", "3"))

saveRDS(seu_int, file = "seu_integrated.rds")
