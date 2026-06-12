import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input)

fig, ax = plt.subplots(figsize=(8, 8))

background = df[~df["is_tumour_high"] & ~df["is_immune"] & ~df["is_interface"]]
tumour = df[df["is_tumour_high"]]
immune = df[df["is_immune"]]
interface = df[df["is_interface"]]

ax.scatter(
    background["x"],
    background["y"],
    s=8,
    c="lightgrey",
    alpha=0.45,
    linewidths=0,
    label="Other spots",
)

ax.scatter(
    tumour["x"],
    tumour["y"],
    s=18,
    c="#D73027",
    alpha=0.75,
    linewidths=0,
    label="Melanoma-high spots",
)

ax.scatter(
    immune["x"],
    immune["y"],
    s=18,
    c="#1A9850",
    alpha=0.75,
    linewidths=0,
    label="Immune spots",
)

ax.scatter(
    interface["x"],
    interface["y"],
    s=26,
    c="#FEE08B",
    alpha=0.95,
    edgecolors="black",
    linewidths=0.2,
    label="Tumour-immune interface",
)

ax.set_title("Tumour-immune interface map")
ax.set_xlabel("Spatial x")
ax.set_ylabel("Spatial y")
ax.invert_yaxis()
ax.set_aspect("equal")
ax.legend(
    frameon=True,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
)

plt.tight_layout()

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(args.output, dpi=300, bbox_inches="tight")
plt.close()
