import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--visium_path", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

visium_path = Path(args.visium_path)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

images = {
    "tissue_hires": visium_path / "spatial" / "tissue_hires_image.png",
    "tissue_lowres": visium_path / "spatial" / "tissue_lowres_image.png",
    "aligned_tissue": visium_path / "spatial" / "aligned_tissue_image.jpg",
    "detected_tissue": visium_path / "spatial" / "detected_tissue_image.jpg",
}

for name, path in images.items():
    if not path.exists():
        continue

    img = mpimg.imread(path)

    plt.figure(figsize=(8, 8))
    plt.imshow(img)
    plt.axis("off")
    plt.title(name.replace("_", " "), fontsize=16)
    plt.tight_layout()
    plt.savefig(outdir / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close()

    # contrast-enhanced version
    img_float = img.astype(float)
    if img_float.max() > 1:
        img_float = img_float / 255.0

    low, high = np.percentile(img_float, [1, 99])
    enhanced = np.clip((img_float - low) / (high - low), 0, 1)

    plt.figure(figsize=(8, 8))
    plt.imshow(enhanced)
    plt.axis("off")
    plt.title(name.replace("_", " ") + " contrast enhanced", fontsize=16)
    plt.tight_layout()
    plt.savefig(outdir / f"{name}_contrast.png", dpi=300, bbox_inches="tight")
    plt.close()
