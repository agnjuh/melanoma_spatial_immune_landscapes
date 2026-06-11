import argparse
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

STATE_COLORS = {
    "inflamed": "#d73027",
    "excluded": "#fc8d59",
    "desert": "#4575b4",
    "immune_niche": "#1a9850",
    "other": "#bdbdbd",
}

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--visium_path", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
visium_path = Path(args.visium_path)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

image_path = visium_path / "spatial" / "tissue_hires_image.png"
img = mpimg.imread(image_path)

x = adata.obs["pxl_col_in_fullres"]
y = adata.obs["pxl_row_in_fullres"]

# The coordinates are full-resolution coordinates.
# The hires image is scaled relative to full resolution.
scale = adata.uns["spatial"][visium_path.name]["scalefactors"]["tissue_hires_scalef"]

x_scaled = x * scale
y_scaled = y * scale

plt.figure(figsize=(8, 8))
plt.imshow(img)

for state, color in STATE_COLORS.items():
    mask = adata.obs["immune_state"] == state
    plt.scatter(
        x_scaled[mask],
        y_scaled[mask],
        s=14,
        c=color,
        label=state,
        linewidths=0,
        alpha=0.78,
    )

plt.title("Immune states overlaid on tissue image")
plt.axis("off")
plt.legend(
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
    fontsize=9,
)

plt.tight_layout()
plt.savefig(outdir / "immune_states_tissue_overlay.png", dpi=300, bbox_inches="tight")
plt.close()


# Marker-score overlays
scores = [
    "melanoma_score",
    "tcell_score",
    "cytotoxic_score",
    "interferon_score",
    "stromal_score",
]

for score in scores:
    plt.figure(figsize=(8, 8))
    plt.imshow(img)
    sc_plot = plt.scatter(
        x_scaled,
        y_scaled,
        c=adata.obs[score],
        s=12,
        linewidths=0,
        alpha=0.80,
    )
    plt.title(score.replace("_", " ") + " overlaid on tissue image")
    plt.axis("off")
    plt.colorbar(sc_plot, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(outdir / f"{score}_tissue_overlay.png", dpi=300, bbox_inches="tight")
    plt.close()
