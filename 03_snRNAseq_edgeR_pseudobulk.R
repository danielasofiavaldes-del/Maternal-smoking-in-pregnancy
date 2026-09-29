# Pseudobulk edgeR likelihood-ratio test, smokers versus non-smokers.
# Design: ~0+condition+batch+gender
# Input: seu_integrated.rds
# Output: deg_pseudobulk.csv

library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(edgeR)

if (!exists("seu_int_fil")) {
  seu_int_fil <- readRDS("seu_integrated.rds")
}

Idents(seu_int_fil) <- "cell_type"
celltypes <- as.character(unique(Idents(seu_int_fil)))
indices <- c("FBp", "EBp")
celltypes <- celltypes[-which(celltypes %in% indices)]

# Raw RNA counts, not the SCT assay. Sum nuclei of the same cell type and sample.
# FBp and EBp are left out.
DefaultAssay(seu_int_fil) <- "RNA"
counts <- GetAssayData(seu_int_fil, slot = "counts", assay = "RNA")
metadata <- seu_int_fil@meta.data
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
ids <- colData(sce)[, c("cell_type", "sample_id")]
summed <- aggregateAcrossCells(sce, ids)

doDEG <- function(celltypes) {
  results <- data.frame()
  for (i in seq_along(celltypes)) {
    label <- celltypes[i]
    current <- summed[, label == summed$cell_type]
    y <- DGEList(counts(current), samples = colData(current))
    keep <- filterByExpr(y, group = current$condition)
    discarded <- current$ncells < 10
    y <- y[, !discarded]
    y <- y[keep, ]
    y <- calcNormFactors(y)
    # Contrast is smokers minus non-smokers, with batch and fetal sex in the model.
    design <- model.matrix(~0 + condition + batch + gender, y$samples)
    rownames(design) <- colnames(y)
    y <- estimateGLMRobustDisp(y, design, verbose = TRUE)
    fit <- glmFit(y, design, robust = TRUE)
    con <- makeContrasts(conditionsmoker - conditioncntrl, levels = design)
    res <- glmLRT(fit, contrast = con)
    fdr <- topTags(res, n = nrow(res))
    res.m <- as.matrix(fdr[[1]])
    file2 <- data.frame(res.m, cluster = celltypes[i])
    file2$gene <- rownames(file2)
    results <- rbind(file2, results)
  }
  results
}

deg_results <- doDEG(celltypes)
write.csv(deg_results, file = "deg_pseudobulk.csv", row.names = FALSE)
