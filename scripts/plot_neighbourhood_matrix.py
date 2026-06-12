import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--value", default="fraction_edges")
args = parser.parse_args()

df = pd.read_csv(args.input)

regions = sorted(set(df["region_a"]).union(set(df["region_b"])))

mat = pd.DataFrame(
    0.0,
    index=regions,
    columns=regions
)

for _, row in df.iterrows():
    mat.loc[row["region_a"], row["region_b"]] = row[args.value]

fig, ax = plt.subplots(figsize=(8, 7))

im = ax.imshow(mat.values)

ax.set_xticks(range(len(mat.columns)))
ax.set_yticks(range(len(mat.index)))

ax.set_xticklabels(mat.columns, rotation=45, ha="right")
ax.set_yticklabels(mat.index)

for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        value = mat.iloc[i, j]
        if value > 0:
            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(args.value)

ax.set_title("Spatial neighbourhood matrix")

plt.tight_layout()

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(args.output, dpi=300, bbox_inches="tight")
plt.close()
