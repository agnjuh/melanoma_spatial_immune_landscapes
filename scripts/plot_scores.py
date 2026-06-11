import argparse
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt


SCORES = [
    "melanoma_score",
    "tcell_score",
    "cytotoxic_score",
    "interferon_score",
    "stromal_score",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    x = adata.obs["pxl_col_in_fullres"]
    y = adata.obs["pxl_row_in_fullres"]

    for score in SCORES:
        plt.figure(figsize=(6, 6))
        plt.scatter(x, -y, c=adata.obs[score], s=8)
        plt.title(score.replace("_", " "))
        plt.axis("equal")
        plt.axis("off")
        plt.colorbar(label=score)
        plt.tight_layout()
        plt.savefig(outdir / f"{score}.png", dpi=300)
        plt.close()


if __name__ == "__main__":
    main()
