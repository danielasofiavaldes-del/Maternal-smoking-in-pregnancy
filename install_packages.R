# Install package versions listed in package_versions.csv.
# Scripts 01-08: R 4.1.2. Script 09 (CellChat 2.1.2): R 4.4.3.

if (getRversion() < "4.1.2") {
  stop("R 4.1.2 or newer is required. The snRNA-seq stack was run in R 4.1.2.")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran <- c(
  Seurat = "4.1.0",
  sctransform = "0.3.3",
  ggplot2 = "3.3.6",
  tidyr = "1.2.0",
  patchwork = "1.1.1",
  cowplot = "1.1.1",
  viridis = "1.6.4",
  RColorBrewer = "1.1.3",
  pheatmap = "1.0.12",
  RVenn = "1.1.0",
  clustree = "0.5.1",
  countsplit = "1.0.0"
)
for (pkg in names(cran)) {
  remotes::install_version(pkg, version = cran[[pkg]], upgrade = "never", quiet = TRUE)
}

bioc_tarballs <- c(
  edgeR = "https://bioconductor.org/packages/3.14/bioc/src/contrib/edgeR_3.36.0.tar.gz",
  limma = "https://bioconductor.org/packages/3.17/bioc/src/contrib/limma_3.56.2.tar.gz",
  SingleCellExperiment = "https://bioconductor.org/packages/3.17/bioc/src/contrib/SingleCellExperiment_1.22.0.tar.gz",
  scuttle = "https://bioconductor.org/packages/3.14/bioc/src/contrib/scuttle_1.4.0.tar.gz",
  slingshot = "https://bioconductor.org/packages/3.14/bioc/src/contrib/slingshot_2.2.1.tar.gz",
  condiments = "https://bioconductor.org/packages/3.14/bioc/src/contrib/condiments_1.2.0.tar.gz",
  EnhancedVolcano = "https://bioconductor.org/packages/3.17/bioc/src/contrib/EnhancedVolcano_1.18.0.tar.gz",
  ComplexHeatmap = "https://bioconductor.org/packages/3.17/bioc/src/contrib/ComplexHeatmap_2.16.0.tar.gz"
)
for (url in bioc_tarballs) {
  remotes::install_url(url, upgrade = "never", quiet = TRUE)
}

remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never", quiet = TRUE)
remotes::install_github("neurorestore/Augur", upgrade = "never", quiet = TRUE)

if (getRversion() >= "4.4.3") {
  remotes::install_version("NMF", version = "0.28", upgrade = "never", quiet = TRUE)
  remotes::install_github("jinworks/CellChat", upgrade = "never", quiet = TRUE)
} else {
  message("CellChat 2.1.2 was run in R 4.4.3. Skip it in this library.")
}
