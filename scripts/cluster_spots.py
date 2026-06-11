import argparse
from pathlib import Path

import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden")

print(adata.obs["leiden"].value_counts().sort_index())

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
adata.write(args.output)
