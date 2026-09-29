#!/bin/bash
# CellBender 0.2.2 remove-background, one sample at a time.
# Input: <sample>/outs/raw_feature_bc_matrix.h5 from 00_cellranger_count.sh
# Output: <sample>_cellbender.rds
#
# 18,000 droplets, 150 epochs, false-positive rate 0.01, CUDA.
# expected-cells 10000 is the value used in the commands.

set -euo pipefail

version=$(cellbender --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
if [[ "${version}" != "0.2.2" ]]; then
  echo "CellBender 0.2.2 is required. Found: ${version:-none}" >&2
  exit 1
fi

mkdir -p cellbender/filtered_h5

run_one() {
  local sample="$1"
  local outdir="cellbender/${sample}"
  mkdir -p "${outdir}"
  cellbender remove-background \
    --input "${sample}/outs/raw_feature_bc_matrix.h5" \
    --output "${outdir}/${sample}_out.h5" \
    --expected-cells 10000 \
    --cuda \
    --total-droplets-included 18000 \
    --fpr 0.01 \
    --epochs 150
  cp "${outdir}/${sample}_out_filtered.h5" "cellbender/filtered_h5/${sample}_out_filtered.h5"
}

for sample in 009 010 011 012 034 035 036 037 038 039 040 041 042 043; do
  run_one "${sample}"
done

# Same object construction as the original reader: min.cells 3, min.features 200.
Rscript -e "
  suppressPackageStartupMessages(library(Seurat))
  files <- list.files('cellbender/filtered_h5', pattern = '_out_filtered\\\\.h5$', full.names = TRUE)
  for (f in files) {
    sample <- sub('_out_filtered.h5$', '', basename(f))
    seu <- CreateSeuratObject(Read10X_h5(f), min.cells = 3, min.features = 200)
    seu\$sample_id <- sample
    seu <- RenameCells(seu, add.cell.id = sample)
    saveRDS(seu, paste0(sample, '_cellbender.rds'))
  }
"
