import argparse
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--figure", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    figure_path = Path(args.figure)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_visium(input_path)

    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True,
        percent_top=None,
        log1p=False,
    )

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].hist(adata.obs["total_counts"], bins=50)
    axes[0].set_title("Total counts")

    axes[1].hist(adata.obs["n_genes_by_counts"], bins=50)
    axes[1].set_title("Genes per spot")

    axes[2].hist(adata.obs["pct_counts_mt"], bins=50)
    axes[2].set_title("Mitochondrial %")

    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()

    adata.write(output_path)


if __name__ == "__main__":
    main()
