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

sample = list(adata.uns["spatial"].keys())[0]

continuous_features = [
    "melanoma_score",
    "tcell_score",
    "cytotoxic_score",
    "interferon_score",
    "stromal_score",
]

for feature in continuous_features:
    sc.pl.spatial(
        adata,
        img_key="hires",
        color=feature,
        library_id=sample,
        alpha_img=0.25,
        spot_size=1.8,
        cmap="viridis",
        show=False,
    )

    plt.title(feature.replace("_", " "), fontsize=18)
    plt.xlabel("")
    plt.ylabel("")
    plt.savefig(
        outdir / f"{feature}_scanpy_spatial.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()


sc.pl.spatial(
    adata,
    img_key="hires",
    color="immune_state",
    library_id=sample,
    alpha_img=0.18,
    spot_size=2.0,
    palette={
        "inflamed": "#d73027",
        "excluded": "#fc8d59",
        "desert": "#4575b4",
        "immune_niche": "#1a9850",
        "other": "#bdbdbd",
    },
    show=False,
)

plt.title("Immune states", fontsize=18)
plt.xlabel("")
plt.ylabel("")
plt.savefig(
    outdir / "immune_state_scanpy_spatial.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()
