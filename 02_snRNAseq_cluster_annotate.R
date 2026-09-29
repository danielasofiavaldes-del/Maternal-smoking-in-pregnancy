# Clustering (resolution 1.0) and the biological cell-type annotation.
# Input: seu_integrated.rds
# Output: annotated seu_integrated.rds, conserved_markers.csv, all_markers.csv, nuclei_proportions.csv

library(Seurat)
library(tidyverse)
library(metap)

if (!exists("seu_int")) {
  seu_int <- readRDS("seu_integrated.rds")
}

DefaultAssay(seu_int) <- "integrated"
set.seed(999999)
seu_int <- RunPCA(seu_int, npcs = 50)
seu_int <- RunUMAP(seu_int, dims = 1:50, reduction = "pca")
seu_int <- FindNeighbors(seu_int, reduction = "pca", dims = 1:30)
seu_int <- FindClusters(
  seu_int,
  resolution = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0),
  random.seed = 999999,
  algorithm = 2
)

Idents(seu_int) <- "integrated_snn_res.1"
seu_int_fil <- subset(seu_int, idents = c("10", "26", "33", "31", "32", "34"), invert = TRUE)

Idents(seu_int_fil) <- "integrated_snn_res.1"
get_conserved <- function(cluster) {
  FindConservedMarkers(
    seu_int_fil,
    ident.1 = cluster,
    grouping.var = "condition",
    only.pos = TRUE
  ) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
conserved_markers <- map_dfr(levels(seu_int_fil), get_conserved)

all_markers <- FindAllMarkers(
  seu_int_fil,
  test.use = "LR",
  min.pct = 0.25,
  min.diff.pct = 0.25,
  max.cells.per.ident = 5000,
  only.pos = TRUE
)

# Annotation of the resolution-1.0 clusters, in level order, from marker expression.
cell_type <- c(
  "CTB", "CTB", "CTBp", "STB", "CTB", "FB", "HBC", "CTB", "HBC", "CTBp", "EB", "STB",
  "CTB", "STB", "EB", "CTB", "CTB", "CTB", "CCT", "PAMM", "STBjuv", "CTBpf", "EB", "MC",
  "EB", "FB", "VEC", "HBCp", "FB"
)
if (length(cell_type) != nlevels(seu_int_fil)) {
  stop(
    "Resolution 1.0 produced ", nlevels(seu_int_fil),
    " clusters; the annotation has ", length(cell_type), " labels."
  )
}
names(cell_type) <- levels(seu_int_fil)
seu_int_fil <- RenameIdents(seu_int_fil, cell_type)
seu_int_fil <- AddMetaData(seu_int_fil, metadata = Idents(seu_int_fil), col.name = "cell_type")
Idents(seu_int_fil) <- "cell_type"

nuclei_per_sample <- prop.table(table(Idents(seu_int_fil), seu_int_fil$sample_id), margin = 2)

saveRDS(seu_int_fil, file = "seu_integrated.rds")
write.csv(conserved_markers, file = "conserved_markers.csv", row.names = FALSE)
write.csv(all_markers, file = "all_markers.csv", row.names = FALSE)
write.csv(as.data.frame(nuclei_per_sample), file = "nuclei_proportions.csv", row.names = FALSE)
