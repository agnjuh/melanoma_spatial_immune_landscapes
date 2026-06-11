import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


FEATURES = {
    "melanoma_score": {
        "title": "Melanoma programme over tissue morphology",
        "label": "Melanoma signature score",
        "filename": "melanoma_programme_morphology_overlay.png",
    },
    "stromal_score": {
        "title": "Stromal programme over tissue morphology",
        "label": "Stromal / ECM signature score",
        "filename": "stromal_programme_morphology_overlay.png",
    },
    "interferon_score": {
        "title": "Interferon response over tissue morphology",
        "label": "Interferon response score",
        "filename": "interferon_response_morphology_overlay.png",
    },
    "tcell_score": {
        "title": "T-cell infiltration over tissue morphology",
        "label": "T-cell signature score",
        "filename": "tcell_infiltration_morphology_overlay.png",
    },
    "cytotoxic_score": {
        "title": "Cytotoxic lymphocyte activity over tissue morphology",
        "label": "Cytotoxic signature score",
        "filename": "cytotoxic_activity_morphology_overlay.png",
    },
}

STATE_COLORS = {
    "inflamed": "#d73027",
    "excluded": "#fc8d59",
    "desert": "#4575b4",
    "immune_niche": "#1a9850",
    "other": "#bdbdbd",
}

STATE_LABELS = {
    "inflamed": "Immune-inflamed",
    "excluded": "Immune-excluded",
    "desert": "Immune-desert",
    "immune_niche": "Immune niche",
    "other": "Unclassified",
}


def clean_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)


def contrast_enhance(img):
    img = img.astype(float)
    if img.max() > 1:
        img = img / 255.0

    low, high = np.percentile(img, [1, 99])
    img = np.clip((img - low) / (high - low), 0, 1)
    return img


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

    image_path = visium_path / "spatial" / "tissue_hires_image.png"
    img = mpimg.imread(image_path)
    img = contrast_enhance(img)

    scale = adata.uns["spatial"][sample]["scalefactors"]["tissue_hires_scalef"]

    x = adata.obs["pxl_col_in_fullres"].to_numpy() * scale
    y = adata.obs["pxl_row_in_fullres"].to_numpy() * scale

    # Tissue-only reference image
    fig, ax = plt.subplots(figsize=(7.5, 7.5), facecolor="white")
    ax.imshow(img)
    clean_axis(ax)
    ax.set_title("Tissue morphology", fontsize=18, pad=12)
    fig.savefig(
        outdir / "tissue_morphology_reference.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    # Continuous score overlays
    for feature, meta in FEATURES.items():
        values = adata.obs[feature].to_numpy()

        vmin = np.nanpercentile(values, 2)
        vmax = np.nanpercentile(values, 98)

        fig, ax = plt.subplots(figsize=(7.8, 7.5), facecolor="white")
        ax.imshow(img, alpha=1.0)

        scplot = ax.scatter(
            x,
            y,
            c=values,
            s=8,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            alpha=0.45,
            linewidths=0,
            edgecolors="none",
        )

        clean_axis(ax)
        ax.set_title(meta["title"], fontsize=17, pad=12)

        cbar = fig.colorbar(scplot, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label(meta["label"], fontsize=10)
        cbar.ax.tick_params(labelsize=9)

        ax.text(
            0.02,
            0.02,
            "Low signal  \u2192  High signal",
            transform=ax.transAxes,
            fontsize=10,
            color="white",
            bbox=dict(
                facecolor="black",
                alpha=0.55,
                edgecolor="none",
                boxstyle="round,pad=0.35",
            ),
        )

        fig.savefig(
            outdir / meta["filename"],
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )
        plt.close(fig)

    # Immune-state overlay
    fig, ax = plt.subplots(figsize=(8.2, 7.5), facecolor="white")
    ax.imshow(img, alpha=1.0)

    plot_order = ["other", "desert", "immune_niche", "inflamed", "excluded"]

    for state in plot_order:
        if state not in adata.obs["immune_state"].astype(str).unique():
            continue

        mask = adata.obs["immune_state"].astype(str) == state

        if state == "other":
            size = 5
            alpha = 0.08
            lw = 0
            edge = "none"
        elif state == "desert":
            size = 8
            alpha = 0.22
            lw = 0
            edge = "none"
        else:
            size = 12
            alpha = 0.72
            lw = 0.10
            edge = "black"

        ax.scatter(
            x[mask],
            y[mask],
            s=size,
            c=STATE_COLORS[state],
            alpha=alpha,
            linewidths=lw,
            edgecolors=edge,
            label=STATE_LABELS[state],
        )

    clean_axis(ax)
    ax.set_title("Spatial immune states over tissue morphology", fontsize=17, pad=12)

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
        outdir / "immune_states_morphology_overlay.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
