"""Deep visual proteomics preprocessing for STB, CTB, and HBC.

Input: report.pg_matrix.tsv, metadata.txt
Output: Spatial_Proteomics_Summary_STB.csv, _CTB.csv, _HBC.csv
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

METADATA_COLUMNS = [
    "patient_id",
    "Smoke_status",
    "Celltype",
    "Celltype_sub",
    "raw_file_id",
    "date",
    "BATCH",
    "protein_hits",
]


def filter_out_str_from_20230315(adata):
    combined_condition = (adata.obs.Celltype_sub == "str") & (adata.obs.date == 20230315)
    return adata[~combined_condition]


def filter_samples_without_hits(adata):
    sc.pp.filter_cells(adata, min_genes=1, inplace=True, copy=False)
    adata.obs["protein_hits"] = adata.obs["n_genes"]
    adata.obs = adata.obs.drop("n_genes", axis=1)
    return adata


def filter_out_hit_threshold(adata, threshold):
    return adata[adata.obs["protein_hits"] > threshold]


def filter_out_hits_and_exclude_ECM_str(adata, threshold):
    return adata[
        (adata.obs["protein_hits"] > threshold)
        | (adata.obs["Celltype_sub"] == "ECM")
        | (adata.obs["Celltype_sub"] == "str")
    ]


def filter_out_contaminants(adata):
    combined_condition = adata.var["Protein.Ids"].str.startswith("Cont_") | adata.var_names.str.startswith("Cont_")
    return adata[:, ~combined_condition]


def filter_invalid_proteins(adata, grouping, threshold):
    valid_protein_names = []
    for group in adata.obs[grouping].unique():
        adata_group = adata[adata.obs[grouping] == group]
        nan_proportions = pd.DataFrame(np.isnan(adata_group.X).mean(axis=0), columns=["nan_proportion"])
        nan_proportions_T_F = nan_proportions["nan_proportion"] <= (1.0 - threshold)
        valid_proteins_group = adata.var[nan_proportions_T_F.values]
        valid_protein_names.extend(list(set(valid_proteins_group.index)))
    valid_proteins_unique = list(set(valid_protein_names))
    return adata[:, valid_proteins_unique]


# Per-protein Gaussian imputation. Width 0.3, downshift -1.8, as in the methods.
def imputation_gaussian(adata, mean_shift=-1.8, std_dev_shift=0.3, perSample=False):
    adata_copy = adata.copy()
    df = pd.DataFrame(data=adata_copy.X, columns=adata_copy.var.index, index=adata_copy.obs_names)
    if perSample:
        df = df.T
    for col in df.columns:
        col_mean = df[col].mean()
        col_std = df[col].std()
        nan_mask = df[col].isnull()
        num_nans = nan_mask.sum()
        random_values = np.random.randn(num_nans)
        shifted_random_values = (col_mean + (mean_shift * col_std)) + (col_std * std_dev_shift) * random_values
        df.loc[nan_mask, col] = shifted_random_values
    if perSample:
        df = df.T
    adata_copy.X = df.values
    return adata_copy


def write_summary_for_limma(adata_batch, celltype):
    column_names = [
        ", ".join(str(value) for value in row)
        for row in adata_batch.obs[METADATA_COLUMNS].to_numpy()
    ]
    summary = pd.DataFrame(
        np.asarray(adata_batch.X).T,
        index=adata_batch.var["Genes"].astype(str),
        columns=column_names,
    )
    summary.index.name = "Protein"
    summary.to_csv(f"Spatial_Proteomics_Summary_{celltype}.csv")


def preprocess_celltype(main_dataset, celltype):
    adata = main_dataset[main_dataset.obs["Celltype"] == celltype].copy()
    adata.X = np.log2(adata.X)
    adata_filtered = filter_invalid_proteins(adata, "Smoke_status", 0.7)
    adata_imputed = imputation_gaussian(adata_filtered, mean_shift=-1.8, std_dev_shift=0.3)
    adata_batch = adata_imputed.copy()
    # ComBat on measurement date. No extra covariates.
    adata_batch.X = sc.pp.combat(adata_imputed, key="BATCH", inplace=False, covariates=None)
    write_summary_for_limma(adata_batch, celltype)
    return adata_batch


# DIA-NN protein groups. Drop CCT, the 15 March 2023 stroma run, shallow samples,
# Cont_ contaminants, and outlier files 2906 (CTB) and 2911 (HBC).
pg_matrix = pd.read_csv("report.pg_matrix.tsv", sep="\t")
rawdata = pg_matrix.iloc[:, 5:].transpose()
sample_metadata = pd.read_csv("metadata.txt", sep="\t")
sample_metadata.index = sample_metadata["Name"]
sample_metadata = sample_metadata.drop("Name", axis=1)
protein_metadata = pg_matrix.iloc[:, :5]
protein_metadata.index = protein_metadata["Protein.Group"]
protein_metadata = protein_metadata.drop("Protein.Group", axis=1)

main_dataset = ad.AnnData(X=rawdata, obs=sample_metadata, var=protein_metadata)
main_dataset = main_dataset[main_dataset.obs.Celltype != "CCT"]
main_dataset = filter_out_str_from_20230315(main_dataset)
main_dataset = filter_samples_without_hits(main_dataset)
main_dataset = filter_out_hit_threshold(main_dataset, 300)
main_dataset = filter_out_hits_and_exclude_ECM_str(main_dataset, 1000)
main_dataset = filter_out_contaminants(main_dataset)
main_dataset = main_dataset[main_dataset.obs.raw_file_id != 2906]
main_dataset = main_dataset[main_dataset.obs.raw_file_id != 2911]

adata_STB = preprocess_celltype(main_dataset, "STB")
adata_CTB = preprocess_celltype(main_dataset, "CTB")
adata_HBC = preprocess_celltype(main_dataset, "HBC")
