import argparse
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

x = adata.obs["pxl_col_in_fullres"]
y = adata.obs["pxl_row_in_fullres"]

plt.figure(figsize=(8, 7))

for region in adata.obs["region_annotation"].cat.categories:
    mask = adata.obs["region_annotation"] == region
    plt.scatter(
        x[mask],
        -y[mask],
        s=9,
        label=region,
        linewidths=0,
    )

plt.title("Annotated spatial tissue regions in melanoma")
plt.axis("equal")
plt.axis("off")
plt.legend(
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
    fontsize=8,
)

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
plt.tight_layout()
plt.savefig(args.output, dpi=300, bbox_inches="tight")
plt.close()
