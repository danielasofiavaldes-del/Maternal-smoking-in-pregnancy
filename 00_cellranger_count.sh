#!/bin/bash
# Cell Ranger 6.1.2 count. Standard GRCh38 reference (exonic reads).
# Writes <sample>/outs/ and <sample>_CR6_unprocessed.rds.
#
# Expects refdata-gex-GRCh38-2020-A/ and fastqs/<sample>/ in this directory.

set -euo pipefail

version=$(cellranger --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
if [[ "${version}" != "6.1.2" ]]; then
  echo "Cell Ranger 6.1.2 is required. Found: ${version:-none}" >&2
  exit 1
fi

transcriptome="refdata-gex-GRCh38-2020-A"
fastqs_dir="fastqs"
samples=(009 010 011 012 034 035 036 037 038 039 040 041 042 043)

for sample in "${samples[@]}"; do
  cellranger count \
    --id="${sample}" \
    --transcriptome="${transcriptome}" \
    --fastqs="${fastqs_dir}/${sample}" \
    --sample="${sample}"

  Rscript -e "
    suppressPackageStartupMessages(library(Seurat))
    sample <- '${sample}'
    mat <- Read10X_h5(file.path(sample, 'outs', 'filtered_feature_bc_matrix.h5'))
    seu <- CreateSeuratObject(counts = mat)
    seu\$sample_id <- sample
    saveRDS(seu, paste0(sample, '_CR6_unprocessed.rds'))
  "
done
