import os
import scanpy as sc
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import numpy as np

OUTDIR = "/mnt/zhangzheng_group/liuz-52/Test_R/01_Figure1_python"
IN_DIR = os.path.join(OUTDIR, "inferCNV_input")
os.makedirs(IN_DIR, exist_ok=True)

h5ad_file = os.path.join(OUTDIR, "Step10_adata_all_final_celltype.h5ad")

adata = sc.read_h5ad(h5ad_file)

print(adata)
print(adata.obs["celltype"].value_counts())

# 优先使用 raw counts
if "counts" in adata.layers:
    X = adata.layers["counts"]
    print("Using adata.layers['counts']")
elif adata.raw is not None:
    X = adata.raw.X
    print("Using adata.raw.X")
else:
    X = adata.X
    print("Warning: using adata.X")

if not sp.issparse(X):
    X = sp.csr_matrix(X)

# inferCNV 需要 genes x cells
X_gene_cell = X.T.tocsr()

# 过滤低表达基因
keep_genes = np.array((X_gene_cell > 0).sum(axis=1)).flatten() >= 10
X_gene_cell = X_gene_cell[keep_genes, :]

genes = adata.var_names[keep_genes].astype(str)
cells = adata.obs_names.astype(str)

mmwrite(os.path.join(IN_DIR, "inferCNV_counts_matrix.mtx"), X_gene_cell)

pd.Series(genes).to_csv(
    os.path.join(IN_DIR, "inferCNV_genes.txt"),
    sep="\t", index=False, header=False
)

pd.Series(cells).to_csv(
    os.path.join(IN_DIR, "inferCNV_cells.txt"),
    sep="\t", index=False, header=False
)

cell_anno = pd.DataFrame({
    "cell": cells,
    "group": adata.obs["celltype"].astype(str).values
})

cell_anno.to_csv(
    os.path.join(IN_DIR, "inferCNV_cell_annotation.txt"),
    sep="\t", index=False, header=False
)

adata.obs[["celltype"]].to_csv(
    os.path.join(IN_DIR, "Step10_cell_metadata_for_inferCNV.csv")
)

print("inferCNV input files saved to:")
print(IN_DIR)
