import argparse
import json
from pathlib import Path

import pandas as pd
import scanpy as sc
import matplotlib.image as mpimg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_10x_mtx(
        input_path / "filtered_feature_bc_matrix",
        var_names="gene_symbols",
        make_unique=True,
    )

    positions = pd.read_csv(input_path / "spatial" / "tissue_positions.csv")
    positions = positions.set_index("barcode")
    adata.obs = adata.obs.join(positions, how="left")

    adata.obsm["spatial"] = adata.obs[
        ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    ].to_numpy()

    with open(input_path / "spatial" / "scalefactors_json.json") as f:
        scalefactors = json.load(f)

    hires = mpimg.imread(input_path / "spatial" / "tissue_hires_image.png")
    lowres = mpimg.imread(input_path / "spatial" / "tissue_lowres_image.png")

    sample_id = input_path.name

    adata.uns["spatial"] = {
        sample_id: {
            "images": {
                "hires": hires,
                "lowres": lowres,
            },
            "scalefactors": scalefactors,
            "metadata": {
                "source": "10x Genomics Human Melanoma IF Stained FFPE"
            },
        }
    }

    adata.write(output_path)
    print(adata)


if __name__ == "__main__":
    main()
