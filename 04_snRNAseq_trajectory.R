# Trophoblast slingshot trajectories, transition genes, and leaf genes.
# Input: seu_integrated.rds
# Output: sce_umap.rds, progression_ks.csv, progression_permutation.csv,
#         pseudotime_quasibinomial.csv, transition_*.csv, leaf_*.csv

library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(slingshot)
library(condiments)
library(countsplit)
library(edgeR)
library(dplyr)
library(tibble)

if (!exists("seu_int_fil")) {
  seu_int_fil <- readRDS("seu_integrated.rds")
}

# Trophoblast only, then a stricter nucleus filter (at least 600 genes).
# countsplit: training half for slingshot, test half for transition and leaf genes.
Idents(seu_int_fil) <- "cell_type"
seu_tropho <- subset(seu_int_fil, idents = c("CTB", "CTBp", "CTBpf", "CCT", "STB", "STBjuv"))
seu_tropho <- subset(seu_tropho, subset = (nFeature_RNA >= 600 & log10GenesPerUMI > 0.8 & percent_mt < 0.25))

DefaultAssay(seu_tropho) <- "RNA"
cnts <- GetAssayData(seu_tropho, slot = "counts")
set.seed(1991)
cnts_split <- countsplit(cnts, epsilon = 0.5)
mat_train <- cnts_split$train
mat_train@x <- as.numeric(mat_train@x)

DefaultAssay(seu_tropho) <- "integrated"
set.seed(1234)
seu_tropho <- RunPCA(seu_tropho, npcs = 30)
set.seed(1234)
seu_tropho <- RunUMAP(seu_tropho, dims = 1:30, reduction = "pca")
set.seed(1234)
seu_tropho <- FindNeighbors(seu_tropho, reduction = "pca", dims = 1:30)
seu_tropho <- Seurat::FindClusters(
  seu_tropho,
  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
  random.seed = 999999,
  algorithm = 2
)

seu_train <- SetAssayData(seu_tropho, slot = "counts", new.data = mat_train)
sce <- as.SingleCellExperiment(seu_train, assay = "RNA")
imbalance <- imbalance_score(
  reducedDims(sce)$UMAP %>% as_tibble() %>% select(UMAP_1, UMAP_2) %>% as.matrix(),
  conditions = colData(sce)$condition,
  k = 20,
  smooth = 40
)

# Cluster 2 expresses proliferation markers and is the root. One curve for both smoking groups.
set.seed(1991)
sce_umap <- slingshot(sce, colData(sce)$integrated_snn_res.0.3, reducedDim = reducedDims(sce)$UMAP, start.clus = "2")

# Smokers versus non-smokers on each lineage, weighted by the slingshot curve weight.
# Kolmogorov-Smirnov, and a permutation of condition labels (10,000 times).
lineage_names <- c("1" = "STB", "2" = "CTB", "3" = "CCT")
progression_ks <- progressionTest(
  sce_umap, conditions = "condition", global = FALSE, lineages = TRUE, method = "KS"
)
progression_permutation <- progressionTest(
  sce_umap, conditions = "condition", global = FALSE, lineages = TRUE,
  method = "Permutation", rep = 10000
)
progression_ks$lineage <- lineage_names[as.character(progression_ks$lineage)]
progression_permutation$lineage <- lineage_names[as.character(progression_permutation$lineage)]

# Quasibinomial model of pseudotime, scaled to (0, 1), on smoking status.
# Same curve weights.
pst <- slingPseudotime(sce_umap)
curve_w <- slingCurveWeights(sce_umap)
condition <- factor(colData(sce_umap)$condition)
quasibinomial_fit <- lapply(seq_len(ncol(pst)), function(i) {
  pt <- pst[, i]
  keep <- !is.na(pt) & curve_w[, i] > 0
  pt <- pt[keep]
  pt01 <- (pt - min(pt)) / (max(pt) - min(pt))
  pt01 <- pmin(pmax(pt01, 1e-6), 1 - 1e-6)
  fit <- glm(
    pseudotime ~ condition,
    family = quasibinomial(),
    weights = curve_w[keep, i],
    data = data.frame(pseudotime = pt01, condition = condition[keep])
  )
  coefs <- as.data.frame(summary(fit)$coefficients)
  coefs$term <- rownames(coefs)
  coefs$lineage <- unname(lineage_names[as.character(i)])
  coefs
})
quasibinomial_fit <- do.call(rbind, quasibinomial_fit)

set.seed(1991)
mat_test <- countsplit(cnts, epsilon = 0.5)$test
mat_test@x <- as.numeric(mat_test@x)
seu_test <- SetAssayData(seu_tropho, slot = "counts", new.data = mat_test)
Idents(seu_test) <- "condition"
seu_cntrl <- subset(seu_test, idents = "cntrl")
seu_cntrl <- FindVariableFeatures(seu_cntrl, selection.method = "vst", nfeatures = 8000)
sce_test <- as.SingleCellExperiment(seu_cntrl, assay = "RNA")

scale_rows <- function(cnts) {
  scaled_cnts <- cnts
  for (i in seq_len(nrow(cnts))) {
    row_indices <- which(cnts[i, ] != 0)
    row_values <- cnts[i, row_indices]
    scaled_cnts[i, row_indices] <- (row_values - min(row_values)) / (max(row_values) - min(row_values))
  }
  scaled_cnts
}

spearman_correlation <- function(cnts, ptime) {
  results_df <- data.frame(gene = rownames(cnts), spearman_corr = NA_real_, p_value = NA_real_, statistic = NA_real_)
  for (i in seq_len(nrow(cnts))) {
    spearman_results <- suppressWarnings(cor.test(cnts[i, ], ptime, method = "spearman"))
    results_df[i, "spearman_corr"] <- spearman_results$estimate
    results_df[i, "p_value"] <- spearman_results$p.value
    results_df[i, "statistic"] <- spearman_results$statistic
  }
  results_df
}

lineage_pseudotime <- function(lineage) {
  cnts_use <- logcounts(sce_test)
  cnts_use <- cnts_use[rownames(cnts_use) %in% VariableFeatures(seu_cntrl), , drop = FALSE]
  sce_use <- sce_umap[, colnames(sce_umap) %in% colnames(cnts_use)]
  sce_use <- sce_use[rownames(sce_use) %in% VariableFeatures(seu_cntrl), ]
  ptime <- slingPseudotime(sce_use)[, lineage]
  ptime <- ptime[!is.na(ptime)]
  list(ptime = ptime, cnts = cnts_use[, colnames(cnts_use) %in% names(ptime), drop = FALSE])
}

transition_one_lineage <- function(lineage) {
  lin <- lineage_pseudotime(lineage)
  ptime <- lin$ptime
  scaled_cnts <- scale_rows(lin$cnts)
  first_cells <- scaled_cnts[, names(ptime)[ptime <= quantile(ptime, 0.2)], drop = FALSE]
  last_cells <- scaled_cnts[, names(ptime)[ptime >= quantile(ptime, 0.8)], drop = FALSE]
  logFC <- data.frame(
    logFC = sparseMatrixStats::rowMeans2(last_cells) / sparseMatrixStats::rowMeans2(first_cells),
    row.names = rownames(scaled_cnts)
  )
  selected_genes <- dplyr::filter(logFC, abs(logFC) > 0.25)
  cnts_selected <- as.matrix(scaled_cnts[rownames(scaled_cnts) %in% rownames(selected_genes), , drop = FALSE])
  corr <- spearman_correlation(cnts_selected, ptime[colnames(cnts_selected)])
  dplyr::filter(corr, abs(spearman_corr) > 0.4)
}

leaf_one_lineage <- function(lineage) {
  lin <- lineage_pseudotime(lineage)
  ptime <- lin$ptime
  sce_leaf <- sce_test[rownames(sce_test) %in% VariableFeatures(seu_cntrl), ]
  colData(sce_leaf)$traj_DEG <- "unknown"
  colData(sce_leaf)$traj_DEG[match(names(ptime)[ptime <= quantile(ptime, 0.2)], colnames(sce_leaf))] <- "start"
  colData(sce_leaf)$traj_DEG[match(names(ptime)[ptime >= quantile(ptime, 0.8)], colnames(sce_leaf))] <- "end"
  sce_leaf$traj_DEG <- factor(sce_leaf$traj_DEG)
  sce_leaf$batch <- factor(sce_leaf$batch)
  sce_leaf$gender <- factor(sce_leaf$gender)
  sce_leaf <- sce_leaf[, sce_leaf$traj_DEG %in% c("start", "end")]
  colData(sce_leaf)$traj_DEG <- droplevels(colData(sce_leaf)$traj_DEG)
  summed <- aggregateAcrossCells(sce_leaf, colData(sce_leaf)[, c("traj_DEG", "sample_id")])
  y <- DGEList(counts(summed), samples = colData(summed))
  y <- calcNormFactors(y)
  design <- model.matrix(~0 + traj_DEG + batch + gender, y$samples)
  y <- estimateGLMRobustDisp(y, design, verbose = TRUE)
  fit <- glmFit(y, design, robust = TRUE)
  con <- makeContrasts(traj_DEGend - traj_DEGstart, levels = design)
  res <- glmLRT(fit, contrast = con)
  as.data.frame(topTags(res, n = nrow(res)))
}

transition_genes <- lapply(1:3, transition_one_lineage)
names(transition_genes) <- c("STB", "CTB", "CCT")
leaf_genes <- lapply(1:3, leaf_one_lineage)
names(leaf_genes) <- c("STB", "CTB", "CCT")

saveRDS(sce_umap, file = "sce_umap.rds")
write.csv(progression_ks, file = "progression_ks.csv", row.names = FALSE)
write.csv(progression_permutation, file = "progression_permutation.csv", row.names = FALSE)
write.csv(quasibinomial_fit, file = "pseudotime_quasibinomial.csv", row.names = FALSE)
for (lineage in names(transition_genes)) {
  write.csv(transition_genes[[lineage]], file = paste0("transition_", lineage, ".csv"), row.names = FALSE)
  write.csv(leaf_genes[[lineage]], file = paste0("leaf_", lineage, ".csv"), row.names = FALSE)
}
