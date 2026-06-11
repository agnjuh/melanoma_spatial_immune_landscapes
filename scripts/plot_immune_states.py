import argparse
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt

COLORS = {
    "inflamed": "#d73027",
    "excluded": "#fc8d59",
    "desert": "#4575b4",
    "immune_niche": "#1a9850",
    "other": "#bdbdbd",
}

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

x = adata.obs["pxl_col_in_fullres"]
y = adata.obs["pxl_row_in_fullres"]

plt.figure(figsize=(7, 7))

for state, color in COLORS.items():
    mask = adata.obs["immune_state"] == state
    plt.scatter(
        x[mask],
        -y[mask],
        s=9,
        c=color,
        label=state,
        linewidths=0,
    )

plt.title("Spatial immune states in melanoma")
plt.axis("equal")
plt.axis("off")
plt.legend(
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
)

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
plt.tight_layout()
plt.savefig(args.output, dpi=300, bbox_inches="tight")
plt.close()
