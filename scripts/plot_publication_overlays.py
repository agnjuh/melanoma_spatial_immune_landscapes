import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


STATE_COLORS = {
    "Immune-inflamed": "#d73027",
    "Immune-excluded": "#fc8d59",
    "Immune-desert": "#4575b4",
    "Immune niche": "#1a9850",
    "Unclassified": "#bdbdbd",
}

STATE_RENAME = {
    "inflamed": "Immune-inflamed",
    "excluded": "Immune-excluded",
    "desert": "Immune-desert",
    "immune_niche": "Immune niche",
    "other": "Unclassified",
}

FEATURES = {
    "melanoma_score": {
        "title": "Melanoma transcriptional programme",
        "label": "Low melanoma signal     High melanoma signal",
    },
    "tcell_score": {
        "title": "T-cell infiltration signature",
        "label": "Low T-cell signal     High T-cell signal",
    },
    "cytotoxic_score": {
        "title": "Cytotoxic lymphocyte activity",
        "label": "Low cytotoxic activity     High cytotoxic activity",
    },
    "interferon_score": {
        "title": "Interferon response",
        "label": "Low IFN response     High IFN response",
    },
    "stromal_score": {
        "title": "Stromal / extracellular matrix programme",
        "label": "Low stromal signal     High stromal signal",
    },
}


def clean_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--visium_path", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    visium_path = Path(args.visium_path)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sample = list(adata.uns["spatial"].keys())[0]

    img = mpimg.imread(visium_path / "spatial" / "tissue_hires_image.png")
    scale = adata.uns["spatial"][sample]["scalefactors"]["tissue_hires_scalef"]

    x = adata.obs["pxl_col_in_fullres"].to_numpy() * scale
    y = adata.obs["pxl_row_in_fullres"].to_numpy() * scale

    background = img.astype(float)
    if background.max() > 1:
        background = background / 255.0

    background = np.clip(background * 0.65, 0, 1)

    for feature, meta in FEATURES.items():
        values = adata.obs[feature].to_numpy()

        vmin = np.nanpercentile(values, 2)
        vmax = np.nanpercentile(values, 98)

        fig, ax = plt.subplots(figsize=(7.5, 7.5), facecolor="white")
        ax.imshow(background)

        scplot = ax.scatter(
            x,
            y,
            c=values,
            s=15,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            alpha=0.78,
            linewidths=0.12,
            edgecolors="black",
        )

        clean_axis(ax)
        ax.set_title(meta["title"], fontsize=17, pad=12)

        ax.text(
            0.02,
            0.02,
            meta["label"],
            transform=ax.transAxes,
            fontsize=9,
            color="white",
            bbox=dict(
                facecolor="black",
                alpha=0.55,
                edgecolor="none",
                boxstyle="round,pad=0.35",
            ),
        )

        cbar = fig.colorbar(scplot, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label("Signature score", fontsize=10)
        cbar.ax.tick_params(labelsize=9)

        fig.savefig(
            outdir / f"{feature}_publication_overlay.png",
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )
        plt.close(fig)

    # Immune state plot
    adata.obs["immune_state_label"] = (
        adata.obs["immune_state"]
        .astype(str)
        .map(STATE_RENAME)
        .fillna("Unclassified")
    )

    fig, ax = plt.subplots(figsize=(8, 7.5), facecolor="white")
    ax.imshow(background)

    plot_order = [
        "Unclassified",
        "Immune-desert",
        "Immune niche",
        "Immune-inflamed",
        "Immune-excluded",
    ]

    for state in plot_order:
        mask = adata.obs["immune_state_label"] == state

        if state == "Unclassified":
            alpha = 0.35
            size = 10
            edge = "none"
            lw = 0
        else:
            alpha = 0.90
            size = 18
            edge = "black"
            lw = 0.15

        ax.scatter(
            x[mask],
            y[mask],
            s=size,
            c=STATE_COLORS[state],
            alpha=alpha,
            linewidths=lw,
            edgecolors=edge,
            label=state,
        )

    clean_axis(ax)
    ax.set_title("Spatial immune states", fontsize=17, pad=12)

    ax.text(
        0.02,
        0.02,
        "Red: immune-inflamed   Blue: immune-desert   Green: immune niche",
        transform=ax.transAxes,
        fontsize=9,
        color="white",
        bbox=dict(
            facecolor="black",
            alpha=0.55,
            edgecolor="none",
            boxstyle="round,pad=0.35",
        ),
    )

    handles, labels = ax.get_legend_handles_labels()
    desired_order = [
        "Immune-inflamed",
        "Immune-excluded",
        "Immune-desert",
        "Immune niche",
        "Unclassified",
    ]

    ordered_handles = []
    ordered_labels = []

    for label in desired_order:
        if label in labels:
            idx = labels.index(label)
            ordered_handles.append(handles[idx])
            ordered_labels.append(labels[idx])

    ax.legend(
        ordered_handles,
        ordered_labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=10,
        markerscale=1.4,
        title="Immune state",
        title_fontsize=11,
    )

    fig.savefig(
        outdir / "immune_states_publication_overlay.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
