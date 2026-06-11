import argparse
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

x = adata.obs["pxl_col_in_fullres"]
y = adata.obs["pxl_row_in_fullres"]

plt.figure(figsize=(7, 7))
for cluster in sorted(adata.obs["leiden"].cat.categories, key=lambda x: int(x)):
    mask = adata.obs["leiden"] == cluster
    plt.scatter(x[mask], -y[mask], s=9, label=cluster, linewidths=0)

plt.title("Spatial Leiden clusters")
plt.axis("equal")
plt.axis("off")
plt.legend(
    title="Cluster",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
)
plt.tight_layout()
plt.savefig(outdir / "spatial_leiden_clusters.png", dpi=300, bbox_inches="tight")
plt.close()

sc.pl.umap(
    adata,
    color=["leiden", "immune_state", "melanoma_score", "tcell_score", "stromal_score"],
    wspace=0.4,
    show=False,
)
plt.savefig(outdir / "umap_clusters_scores.png", dpi=300, bbox_inches="tight")
plt.close()
