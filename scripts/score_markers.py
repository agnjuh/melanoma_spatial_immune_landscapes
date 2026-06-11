import argparse
from pathlib import Path

import scanpy as sc
import pandas as pd


MARKER_SETS = {
    "melanoma_score": ["MLANA", "PMEL", "TYR", "MITF", "SOX10"],
    "tcell_score": ["CD3D", "CD3E", "CD8A", "TRAC", "TRBC1"],
    "cytotoxic_score": ["NKG7", "GZMB", "PRF1", "GNLY"],
    "interferon_score": ["IFNG", "CXCL9", "CXCL10", "STAT1", "IRF1"],
    "stromal_score": ["COL1A1", "COL3A1", "FN1", "FAP", "TGFB1"],
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    for score_name, genes in MARKER_SETS.items():
        present = [g for g in genes if g in adata.var_names]
        print(score_name, present)
        if present:
            sc.tl.score_genes(adata, gene_list=present, score_name=score_name)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write(args.output)
    print(adata)


if __name__ == "__main__":
    main()
