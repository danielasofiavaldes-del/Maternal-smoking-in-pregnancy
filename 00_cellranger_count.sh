#!/bin/bash
# Cell Ranger 6.1.2 count.
# Reference: GRCh38 pre-mRNA 3.0.0 (transcript features reannotated as exons).
# Writes <sample>/outs/ and <sample>_CR6_unprocessed.rds.
#
# Expects GRCh38_premRNA_3.0.0/ and fastqs/<sample>/ in this directory.

set -euo pipefail

version=$(cellranger --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
if [[ "${version}" != "6.1.2" ]]; then
  echo "Cell Ranger 6.1.2 is required. Found: ${version:-none}" >&2
  exit 1
fi

transcriptome="GRCh38_premRNA_3.0.0"
fastqs_dir="fastqs"
samples=(009 010 011 012 034 035 036 037 038 039 040 041 042 043)

# Samples 038, 040 and 043 are counted here and dropped later, during QC.
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
