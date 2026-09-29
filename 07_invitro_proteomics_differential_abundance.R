# In-vitro STB limma (CSE versus control at day 6 and day 7) and the untreated differentiation correlation.
# Input: Bulk_Proteomics_Summary.csv
# Output: da_stb_d6.csv, da_stb_d7.csv, differentiation_correlation.csv

library(SingleCellExperiment)
library(limma)

protein_data <- read.csv("Bulk_Proteomics_Summary.csv", header = TRUE)
sample_names <- colnames(protein_data)[2:ncol(protein_data)]
metadata <- data.frame(matrix(unlist(strsplit(sample_names, "_")), nrow = length(sample_names), byrow = TRUE))
colnames(metadata) <- c("CellType", "Timepoint", "Treatment", "Passage", "ID")
counts <- as.matrix(protein_data[, 2:ncol(protein_data)])
rownames(counts) <- protein_data$Protein
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

# Column names are CellType_Timepoint_Treatment_Passage_ID.
# d524h is day 6 and d548h is day 7. Passage is in the design.
limma_cse <- function(time_point) {
  stb_sce <- sce[, colData(sce)$CellType == "STB" & colData(sce)$Timepoint == time_point]
  design <- model.matrix(~ 0 + Treatment + Passage, data = colData(stb_sce))
  fit <- lmFit(counts(stb_sce), design)
  contrast_matrix <- makeContrasts(TreatmentCSE - Treatmentcontrol, levels = design)
  fit_contrast <- contrasts.fit(fit, contrast_matrix)
  fit_ebayes <- eBayes(fit_contrast)
  results <- topTable(fit_ebayes, adjust.method = "BH", number = Inf)
  results$Protein <- rownames(results)
  results
}

da_stb_d6 <- limma_cse("d524h")
da_stb_d7 <- limma_cse("d548h")

control_metadata <- metadata[metadata$Treatment == "control", ]
control_counts <- counts[, metadata$Treatment == "control"]
diff_order <- c("d0", "d5", "d524h", "d548h")
timepoint_mapping <- match(control_metadata$Timepoint, diff_order)
ordered_control_counts <- control_counts[, order(timepoint_mapping)]
differentiation_correlation <- apply(ordered_control_counts, 1, function(x) cor(x, timepoint_mapping))
differentiation_proteins <- names(differentiation_correlation[abs(differentiation_correlation) > 0.5])

write.csv(da_stb_d6, file = "da_stb_d6.csv", row.names = FALSE)
write.csv(da_stb_d7, file = "da_stb_d7.csv", row.names = FALSE)
write.csv(
  data.frame(
    Protein = differentiation_proteins,
    correlation = differentiation_correlation[differentiation_proteins]
  ),
  file = "differentiation_correlation.csv",
  row.names = FALSE
)
