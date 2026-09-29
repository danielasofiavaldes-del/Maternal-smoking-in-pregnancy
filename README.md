# Maternal smoking in early pregnancy

Code for the computational analyses in:

Valdes DS, Nimo J, Nonn O, Ulrich J, Forstner D, Haider S, Ressler M, Pignitter M, Obermayer-Pietsch B, Müller DN, Dechend R, Gauster M, Coscia F, Herse F. Maternal smoking in early pregnancy disrupts placental function through syncytiotrophoblast and macrophage dysregulation. bioRxiv 2025.07.10.664172. <https://doi.org/10.1101/2025.07.10.664172>

First-trimester placental villi from cotinine-confirmed smokers and non-smokers were profiled by single-nucleus RNA-seq and by deep visual proteomics of the same tissues. Selected changes were checked in an independent cohort and in primary trophoblast cells exposed to cigarette-smoke extract. The syncytiotrophoblast and Hofbauer cells carry most of the smoking-associated signal.

## Data

- snRNA-seq: [GSE304028](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304028)
- Proteomics: [PXD066340](https://www.ebi.ac.uk/pride/archive/projects/PXD066340)

## Code

Each scrip included herein reads the file written by the previous step. `00_cellranger_count.sh` needs Cell Ranger 6.1.2. `00b_cellbender_remove_background.sh` needs CellBender 0.2.2 and a GPU (`--cuda`). The other R scripts were run in R 4.1.2, except CellChat (`09`), which was run in R 4.4.3. Versions are listed in `package_versions.csv`.

The scripts cover counting against the GRCh38 pre-mRNA 3.0.0 reference, CellBender ambient correction (18,000 droplets, 150 epochs, false-positive rate 0.01), nucleus QC, doublet scoring (the calls are kept, the nuclei stay in the object), SCTransform integration without regressing mitochondrial counts or sequencing depth, clustering and cell-type annotation, pseudobulk edgeR, the trophoblast slingshot trajectory, proteomics preprocessing and limma, CellChat, and the computational figure panels. Short notes on the parameters are in the scripts.

Perseus enrichment, Cytoscape networks, qPCR, ELISA, imaging, and the mitochondrial toxicity assay were done in other software.
