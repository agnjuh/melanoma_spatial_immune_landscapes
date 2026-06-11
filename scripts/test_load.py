import scanpy as sc

adata = sc.read_10x_mtx(
    "data/raw/melanoma_if_ffpe/filtered_feature_bc_matrix",
    var_names="gene_symbols",
    make_unique=True
)

print(adata)
