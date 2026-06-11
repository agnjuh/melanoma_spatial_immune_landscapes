import argparse
from pathlib import Path

import scanpy as sc

ANNOTATIONS = {
    "0": "melanoma_rich_tumour",
    "1": "collagen_rich_stroma",
    "2": "mixed_tumour_stromal",
    "3": "intermediate_tumour_region",
    "4": "vascular_stroma",
    "5": "basal_epithelial",
    "6": "differentiated_epidermal",
    "7": "fibroblast_matrix_stroma",
    "8": "inflammatory_epidermal",
    "9": "basal_epidermal",
    "10": "differentiated_epidermal_2",
}

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

adata.obs["region_annotation"] = adata.obs["leiden"].map(ANNOTATIONS).astype("category")

print(adata.obs["region_annotation"].value_counts())

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
adata.write(args.output)
