import argparse
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

counts = adata.obs["immune_state"].value_counts().sort_index()

summary = pd.DataFrame({
    "state": counts.index,
    "n_spots": counts.values
})

summary["percent"] = (
    summary["n_spots"] / summary["n_spots"].sum() * 100
).round(2)

summary.to_csv(args.output, index=False)

print(summary)
