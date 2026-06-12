import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

PROGRAMMES = [
    "melanoma_score",
    "tcell_score",
    "cytotoxic_score",
    "interferon_score",
    "stromal_score",
]

LABELS = {
    "melanoma_score": "Melanoma",
    "tcell_score": "T cell",
    "cytotoxic_score": "Cytotoxic",
    "interferon_score": "Interferon",
    "stromal_score": "Stromal",
}

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input)

corr = pd.DataFrame(
    np.eye(len(PROGRAMMES)),
    index=PROGRAMMES,
    columns=PROGRAMMES
)

for _, row in df.iterrows():
    a = row["feature_a"]
    b = row["feature_b"]
    r = row["spearman_r"]

    corr.loc[a, b] = r
    corr.loc[b, a] = r

corr.index = [LABELS[x] for x in corr.index]
corr.columns = [LABELS[x] for x in corr.columns]

fig, ax = plt.subplots(figsize=(7, 6))

im = ax.imshow(corr.values, vmin=-1, vmax=1)

ax.set_xticks(range(len(corr.columns)))
ax.set_yticks(range(len(corr.index)))

ax.set_xticklabels(corr.columns, rotation=45, ha="right")
ax.set_yticklabels(corr.index)

for i in range(len(corr)):
    for j in range(len(corr)):
        ax.text(
            j,
            i,
            f"{corr.iloc[i,j]:.2f}",
            ha="center",
            va="center",
            fontsize=9,
        )

cbar = plt.colorbar(im)
cbar.set_label("Spearman correlation")

ax.set_title("Programme correlation matrix")

plt.tight_layout()

Path(args.output).parent.mkdir(parents=True, exist_ok=True)

plt.savefig(
    args.output,
    dpi=300,
    bbox_inches="tight"
)

plt.close()
