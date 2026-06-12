import argparse
from itertools import combinations
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, pearsonr


PROGRAMMES = [
    "melanoma_score",
    "tcell_score",
    "cytotoxic_score",
    "interferon_score",
    "stromal_score",
]


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

rows = []

for feature_a, feature_b in combinations(PROGRAMMES, 2):
    x = adata.obs[feature_a]
    y = adata.obs[feature_b]

    spearman_r, spearman_p = spearmanr(x, y)
    pearson_r, pearson_p = pearsonr(x, y)

    rows.append({
        "feature_a": feature_a,
        "feature_b": feature_b,
        "spearman_r": round(spearman_r, 4),
        "spearman_p": spearman_p,
        "pearson_r": round(pearson_r, 4),
        "pearson_p": pearson_p,
    })

df = pd.DataFrame(rows)

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(args.output, index=False)

print(df)
