import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--cluster_key", default="immune_state")
parser.add_argument("--radius", type=float, default=250.0)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

if args.cluster_key not in adata.obs.columns:
    raise ValueError(f"Column not found in adata.obs: {args.cluster_key}")

if "spatial" not in adata.obsm:
    raise ValueError("Missing adata.obsm['spatial'].")

coords = np.asarray(adata.obsm["spatial"])
labels = adata.obs[args.cluster_key].astype(str).values

tree = cKDTree(coords)
pairs = tree.query_pairs(r=args.radius)

rows = []

for i, j in pairs:
    a = labels[i]
    b = labels[j]

    rows.append({"region_a": a, "region_b": b})
    rows.append({"region_a": b, "region_b": a})

edges = pd.DataFrame(rows)

if edges.empty:
    raise ValueError("No neighbouring pairs detected. Try increasing --radius.")

counts = (
    edges
    .groupby(["region_a", "region_b"])
    .size()
    .reset_index(name="n_edges")
)

total_edges = counts["n_edges"].sum()
counts["fraction_edges"] = counts["n_edges"] / total_edges

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
counts.to_csv(args.output, index=False)

print(counts.sort_values("n_edges", ascending=False).head(20))
