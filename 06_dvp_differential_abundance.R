# Limma differential abundance, smokers versus non-smokers, after averaging by patient.
# Input: Spatial_Proteomics_Summary_STB.csv, _CTB.csv, _HBC.csv
# Output: da_STB.csv, da_CTB.csv, da_HBC.csv

library(SingleCellExperiment)
library(limma)
library(S4Vectors)

limma_smoke <- function(summary_file) {
  protein_data <- read.csv(summary_file, header = TRUE, check.names = FALSE)
  metadata_colnames <- c(
    "patient_id", "Smoke_status", "Celltype", "Celltype_sub",
    "raw_file_id", "date", "BATCH", "protein_hits"
  )
  long_names <- colnames(protein_data)[2:length(protein_data)]
  metadata_df <- as.data.frame(do.call(rbind, strsplit(long_names, ", ")))
  colnames(metadata_df) <- metadata_colnames
  sample_names <- sub("^(?:[^,]+,){4}\\s*(\\d+).*", "\\1", long_names)

  counts <- as.matrix(protein_data[2:length(protein_data)])
  colnames(counts) <- sample_names
  rownames(counts) <- protein_data$Protein
  sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata_df)

  patient_ids <- colData(sce)$patient_id
  group_indices <- split(seq_len(ncol(counts(sce))), patient_ids)
  averaged_counts <- sapply(group_indices, function(indices) {
    rowMeans(counts(sce)[, indices, drop = FALSE])
  })

  patient_meta <- metadata_df[!duplicated(metadata_df$patient_id), ]
  patient_meta <- patient_meta[match(colnames(averaged_counts), patient_meta$patient_id), ]
  sce_avg <- SingleCellExperiment(
    assays = list(counts = averaged_counts),
    colData = DataFrame(
      patient_id = patient_meta$patient_id,
      Smoke_status = patient_meta$Smoke_status
    )
  )

  # Metadata labels are Smoke and No_smoke. Replicates are already averaged by patient.
  design <- model.matrix(~ 0 + Smoke_status, data = colData(sce_avg))
  fit <- lmFit(counts(sce_avg), design)
  contrast_matrix <- makeContrasts(Smoke_statusSmoke - Smoke_statusNo_smoke, levels = design)
  fit_contrast <- contrasts.fit(fit, contrast_matrix)
  fit_ebayes <- eBayes(fit_contrast)
  results <- topTable(fit_ebayes, adjust.method = "BH", number = Inf)
  results$Protein <- rownames(results)
  results
}

da_STB <- limma_smoke("Spatial_Proteomics_Summary_STB.csv")
da_CTB <- limma_smoke("Spatial_Proteomics_Summary_CTB.csv")
da_HBC <- limma_smoke("Spatial_Proteomics_Summary_HBC.csv")

write.csv(da_STB, file = "da_STB.csv", row.names = FALSE)
write.csv(da_CTB, file = "da_CTB.csv", row.names = FALSE)
write.csv(da_HBC, file = "da_HBC.csv", row.names = FALSE)
