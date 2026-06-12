import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input)

immune = df[df["is_immune"]].copy()

if immune.empty:
    raise ValueError("No immune spots found in interface distance table.")

fig, ax = plt.subplots(figsize=(7, 5))

ax.hist(
    immune["distance_to_nearest_tumour"],
    bins=40,
    edgecolor="black",
    linewidth=0.4,
)

median_distance = immune["distance_to_nearest_tumour"].median()
ax.set_xlabel("Distance to nearest melanoma-high spot (spatial coordinate units)")
ax.axvline(
    median_distance,
    linestyle="--",
    linewidth=2,
)

ax.set_title("Distance from immune spots to nearest melanoma-high region")
ax.set_xlabel("Distance to nearest melanoma-high spot (spatial coordinate units)")
ax.set_ylabel("Number of immune spots")

ax.text(
    0.98,
    0.95,
    f"Median = {median_distance:.1f}",
    transform=ax.transAxes,
    ha="right",
    va="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
)

plt.tight_layout()

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(args.output, dpi=300, bbox_inches="tight")
plt.close()
