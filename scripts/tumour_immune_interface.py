import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--summary", required=True)
parser.add_argument("--melanoma_quantile", type=float, default=0.70)
parser.add_argument("--interface_distance", type=float, default=150.0)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

required_obs = [
    "melanoma_score",
    "immune_state",
]

for col in required_obs:
    if col not in adata.obs.columns:
        raise ValueError(f"Missing required adata.obs column: {col}")

if "spatial" not in adata.obsm:
    raise ValueError("Missing adata.obsm['spatial'] coordinates.")

coords = np.asarray(adata.obsm["spatial"])

melanoma_threshold = adata.obs["melanoma_score"].quantile(args.melanoma_quantile)

tumour_mask = adata.obs["melanoma_score"] >= melanoma_threshold
immune_mask = adata.obs["immune_state"].isin(["inflamed", "immune_niche"])

tumour_coords = coords[tumour_mask.values]
immune_coords = coords[immune_mask.values]

if tumour_coords.shape[0] == 0:
    raise ValueError("No tumour-high spots detected. Try lowering --melanoma_quantile.")

if immune_coords.shape[0] == 0:
    raise ValueError("No immune spots detected using states: inflamed, immune_niche.")

tumour_tree = cKDTree(tumour_coords)
immune_tree = cKDTree(immune_coords)

distance_to_tumour = tumour_tree.query(coords, k=1)[0]
distance_to_immune = immune_tree.query(coords, k=1)[0]

interface_mask = (
    tumour_mask.values | immune_mask.values
) & (
    distance_to_tumour <= args.interface_distance
) & (
    distance_to_immune <= args.interface_distance
)

out = pd.DataFrame({
    "spot_id": adata.obs_names,
    "x": coords[:, 0],
    "y": coords[:, 1],
    "melanoma_score": adata.obs["melanoma_score"].values,
    "immune_state": adata.obs["immune_state"].astype(str).values,
    "is_tumour_high": tumour_mask.values,
    "is_immune": immune_mask.values,
    "distance_to_nearest_tumour": distance_to_tumour,
    "distance_to_nearest_immune": distance_to_immune,
    "is_interface": interface_mask,
})

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
Path(args.summary).parent.mkdir(parents=True, exist_ok=True)

out.to_csv(args.output, index=False)

summary = pd.DataFrame([
    {
        "metric": "n_spots",
        "value": adata.n_obs,
    },
    {
        "metric": "melanoma_quantile",
        "value": args.melanoma_quantile,
    },
    {
        "metric": "melanoma_score_threshold",
        "value": melanoma_threshold,
    },
    {
        "metric": "n_tumour_high_spots",
        "value": int(tumour_mask.sum()),
    },
    {
        "metric": "n_immune_spots",
        "value": int(immune_mask.sum()),
    },
    {
        "metric": "n_interface_spots",
        "value": int(interface_mask.sum()),
    },
    {
        "metric": "percent_interface_spots",
        "value": round(float(interface_mask.mean() * 100), 2),
    },
    {
        "metric": "median_immune_to_tumour_distance",
        "value": round(float(np.median(distance_to_tumour[immune_mask.values])), 2),
    },
    {
        "metric": "mean_immune_to_tumour_distance",
        "value": round(float(np.mean(distance_to_tumour[immune_mask.values])), 2),
    },
    {
        "metric": "p90_immune_to_tumour_distance",
        "value": round(float(np.percentile(distance_to_tumour[immune_mask.values], 90)), 2),
    },
    {
        "metric": "interface_distance_threshold",
        "value": args.interface_distance,
    },
])

summary.to_csv(args.summary, index=False)

print(summary)
