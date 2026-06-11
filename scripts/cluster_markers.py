import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon"
)

Path(args.outdir).mkdir(parents=True, exist_ok=True)

clusters = adata.obs["leiden"].cat.categories

for cluster in clusters:
    df = sc.get.rank_genes_groups_df(
        adata,
        group=cluster
    )

    df.to_csv(
        f"{args.outdir}/cluster_{cluster}_markers.csv",
        index=False
    )

    print("\n")
    print(f"=== Cluster {cluster} ===")
    print(df.head(10)[["names", "logfoldchanges", "pvals_adj"]])

adata.write(
    f"{args.outdir}/annotated_clusters.h5ad"
)
