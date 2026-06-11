import argparse
import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

tumor_thr = adata.obs["melanoma_score"].quantile(0.60)
tcell_thr = adata.obs["tcell_score"].quantile(0.70)
ifn_thr = adata.obs["interferon_score"].quantile(0.70)
stroma_thr = adata.obs["stromal_score"].quantile(0.70)

state = np.array(["other"] * adata.n_obs, dtype=object)

inflamed = (
    (adata.obs["melanoma_score"] > tumor_thr) &
    ((adata.obs["tcell_score"] > tcell_thr) |
     (adata.obs["interferon_score"] > ifn_thr))
)

excluded = (
    (adata.obs["melanoma_score"] > tumor_thr) &
    (adata.obs["stromal_score"] > stroma_thr) &
    (adata.obs["tcell_score"] <= tcell_thr)
)

desert = (
    (adata.obs["melanoma_score"] > tumor_thr) &
    (adata.obs["tcell_score"] <= tcell_thr) &
    (adata.obs["interferon_score"] <= ifn_thr) &
    (adata.obs["stromal_score"] <= stroma_thr)
)

immune_niche = (
    (adata.obs["tcell_score"] > tcell_thr) &
    (adata.obs["melanoma_score"] <= tumor_thr)
)

state[desert] = "desert"
state[excluded] = "excluded"
state[immune_niche] = "immune_niche"
state[inflamed] = "inflamed"

adata.obs["immune_state"] = pd.Categorical(
    state,
    categories=["inflamed", "excluded", "desert", "immune_niche", "other"]
)

print(adata.obs["immune_state"].value_counts(dropna=False))

adata.write(args.output)
