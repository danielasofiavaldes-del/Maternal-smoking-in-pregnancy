# Nucleus QC, doublet removal, and SCTransform integration.
# Input: *_CR6_unprocessed.rds and cell_cycle.rda
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

read_single_cell_data <- function(input_data_pattern = "_CR6_unprocessed.rds") {
  list <- list.files(pattern = "_CR6_unprocessed\\.rds$", full.names = TRUE)
  list_names <- gsub("SP082_", "", basename(list))
  list_names <- gsub("_SeuratInitial_20210527.rds", "", list_names)
  list_names <- gsub("_CR6_unprocessed.rds", "", list_names)
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

FindDoublets <- function(list) {
  set.seed(1991)
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
    seu <- FindClusters(seu, resolution = 0.5)
    sweep.res.list <- paramSweep_v3(seu, PCs = 1:min.pc, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
    optimal.pk <- bcmvn.max$pK
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    nExp_poi <- round(optimal.pk * nrow(seu@meta.data))
    seu <- doubletFinder_v3(seu, PCs = 1:min.pc, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    DF.name <- colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]
    seu <- seu[, seu@meta.data[, DF.name] == "Singlet"]
    seu
  })
}

integrate_SCT_single_cell_data <- function(seul, cycle_genes) {
  set.seed(1234)
  load(cycle_genes)
  for (i in seq_along(seul)) {
    seul[[i]] <- CellCycleScoring(seul[[i]], g2m.features = g2m_genes, s.features = s_genes)
  }
  features <- SelectIntegrationFeatures(object.list = seul, nfeatures = 4000)
  seul <- PrepSCTIntegration(object.list = seul, anchor.features = features)
  integ_anchors <- FindIntegrationAnchors(object.list = seul, anchor.features = features, normalization.method = "SCT")
  seu_int <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
  DefaultAssay(seu_int) <- "integrated"
  seu_int
}

sample_ids <- c("009", "010", "043", "036", "035", "042", "040", "011", "012", "039", "041", "037", "034", "038")

data_list <- read_single_cell_data()
data_list <- Preprocessing(data_list)
data_list <- filter(data_list)
# Samples 038, 040 and 043, dropped before integration in the CellBender script
data_list <- data_list[!vapply(data_list, function(seu) unique(seu$sample_id) %in% c("038", "040", "043"), logical(1))]
datal_filtered_nd <- FindDoublets(data_list)

seu_int <- integrate_SCT_single_cell_data(datal_filtered_nd, cycle_genes = "cell_cycle.rda")

seu_int$GA <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("68", "70", "75", "49", "52", "61", "65", "74", "76", "35", "43", "58", "60", "64"))
seu_int$condition <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("smoker", "smoker", "smoker", "smoker", "smoker", "smoker", "smoker", "cntrl", "cntrl", "cntrl", "cntrl", "cntrl", "cntrl", "cntrl"))
seu_int$gender <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("female", "female", "male", "female", "male", "male", "female", "female", "female", "female", "female", "male", "male", "female"))
seu_int$batch <- plyr::mapvalues(seu_int$sample_id, from = sample_ids, to = c("1", "1", "4", "3", "2", "4", "3", "1", "1", "3", "3", "3", "2", "3"))

saveRDS(seu_int, file = "seu_integrated.rds")
