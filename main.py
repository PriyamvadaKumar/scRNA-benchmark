#!/usr/bin/env python


import os
import time
import json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

import umap
import pymde
from sklearn.manifold import TSNE

# CONFIG

DATASETS = {
    "pbmc3k": {
        "path": "pbmc3k.h5ad",
        "min_counts": 100,
        "min_genes": 100,
        "cell_type_col": "cell_type",
    },
    "pbmc68k": {
        "path": "pbmc68k.h5ad",
        "min_counts": 1,
        "min_genes": 1,
        "cell_type_col": "celltype",
    },
}

N_PCS = 100
N_NEIGHBORS = 15
TOP_VAR_GENES = 3000
RANDOM_STATE = 42
MIN_CELLS_PER_GENE = 3

os.makedirs("results_pbmc3k", exist_ok=True)
os.makedirs("results_pbmc68k", exist_ok=True)

sns.set_style("white")
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, figsize=(5, 4))


# HELPERS
def ensure_cell_type(adata, cfg):
    col = cfg["cell_type_col"]
    print(f"  Looking for cell type column: '{col}'")
    if col in adata.obs.columns:
        adata.obs["cell_type"] = adata.obs[col].astype("category")
        print(f"  Found '{col}' → {adata.obs['cell_type'].nunique()} cell types")
    else:
        available = adata.obs.columns.tolist()
        raise KeyError(f"Cell type column '{col}' not found! Available: {available}")
    return adata


def is_raw_counts(adata):
    """Auto-detect if data is raw counts or pre-normalized."""
    X_sample = adata.X[:100, :100]
    if hasattr(X_sample, "toarray"):
        X_sample = X_sample.toarray()
    values = X_sample.flatten()
    values = values[values > 0]

    if len(values) == 0:
        return True  # fallback

    is_integer = np.issubdtype(X_sample.dtype, np.integer) or np.all(values == values.astype(int))
    high_max = values.max() > 100

    return is_integer and high_max


def run_embedding_full(X_pca, method, n_cells):
    start = time.time()
    if method == "umap":
        n = min(N_NEIGHBORS, n_cells - 1)
        reducer = umap.UMAP(n_neighbors=n, min_dist=0.1, n_components=2, random_state=RANDOM_STATE)
        emb = reducer.fit_transform(X_pca)
    elif method == "mde":
        n = min(N_NEIGHBORS, n_cells - 1)
        mde = pymde.preserve_neighbors(X_pca, embedding_dim=2, n_neighbors=n, device="cpu")
        emb = mde.embed().cpu().numpy()
    elif method == "tsne":
        perplexity = min(30, n_cells // 5)
        tsne = TSNE(n_components=2, perplexity=perplexity, random_state=RANDOM_STATE, n_jobs=-1, init="random")
        emb = tsne.fit_transform(X_pca)
    return emb, time.time() - start


def clean_and_preprocess(adata_raw, cfg):
    """Robust preprocessing with auto-detection of raw vs normalized."""
    
    # 1. Filter genes
    sc.pp.filter_genes(adata_raw, min_cells=MIN_CELLS_PER_GENE)
    total_counts = np.ravel(adata_raw.X.sum(axis=1))
    keep = total_counts > 0
    if not keep.any():
        raise ValueError("All cells have zero counts after QC!")
    adata_raw = adata_raw[keep].copy()
    print(f"  After initial QC: {adata_raw.n_obs} cells, {adata_raw.n_vars} genes")

    # 2. HVG on current data
    n_top = min(TOP_VAR_GENES, adata_raw.n_vars // 2)
    print(f"  HVG: Seurat v1 (top {n_top})")
    sc.pp.highly_variable_genes(
        adata_raw,
        n_top_genes=n_top,
        flavor="seurat",
        subset=True,
        inplace=True,
    )
    print(f"  → {adata_raw.n_vars} HVGs selected")

    # 3. Remove zero-count cells AFTER HVG
    total_counts_hvg = np.ravel(adata_raw.X.sum(axis=1))
    keep = total_counts_hvg > 0
    if not keep.any():
        raise ValueError("All cells have zero counts after HVG subset!")
    adata = adata_raw[keep].copy()
    print(f"  After HVG zero-count removal: {adata.n_obs} cells")

    # 4. AUTO-DETECT: raw or pre-normalized?
    is_raw = is_raw_counts(adata)
    print(f"  DATA TYPE: {'RAW COUNTS' if is_raw else 'PRE-NORMALIZED'}")

    # 5. Normalize only if raw
    if is_raw:
        print("  Applying normalization + log1p")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("  Skipping normalization (already normalized)")

    # 6. Final safety: remove NaN/Inf
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    mask = np.isfinite(X).all(axis=1)
    if not mask.all():
        n_remove = (~mask).sum()
        print(f"  Removing {n_remove} cells with NaN/Inf")
        adata = adata[mask].copy()

    if adata.n_obs == 0:
        raise ValueError("No cells remaining after preprocessing!")

    # 7. PCA
    n_pcs = min(N_PCS, adata.n_obs - 1, adata.n_vars - 1)
    if n_pcs <= 0:
        raise ValueError(f"Invalid n_pcs={n_pcs}")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=RANDOM_STATE)
    print(f"  PCA: {n_pcs} components")
    return adata


# MAIN LOOP

all_results = []

for ds_name, cfg in DATASETS.items():
    print(f"\n{'='*60}")
    print(f"RUNNING {ds_name.upper()}")
    print(f"{'='*60}")

    OUT_DIR = f"results_{ds_name}"
    os.makedirs(OUT_DIR, exist_ok=True)

    adata_raw = sc.read_h5ad(cfg["path"])
    print(f"Loaded: {adata_raw.n_obs:,} cells, {adata_raw.n_vars:,} genes")

    # Custom QC for PBMC3k
    if ds_name == "pbmc3k":
        sc.pp.filter_cells(adata_raw, min_counts=cfg["min_counts"])
        sc.pp.filter_cells(adata_raw, min_genes=cfg["min_genes"])

    # Cell type
    adata_raw = ensure_cell_type(adata_raw, cfg)

    try:
        adata = clean_and_preprocess(adata_raw, cfg)

        true_labels = adata.obs["cell_type"].cat.codes.values
        X_pca = adata.obsm["X_pca"]

        results = {
            "dataset": ds_name,
            "n_cells": adata.n_obs,
            "n_types": adata.obs["cell_type"].nunique(),
            "methods": {},
        }
        embeddings = {}

        for method in ["umap", "mde", "tsne"]:
            print(f"  → {method.upper()} ...")
            emb, runtime = run_embedding_full(X_pca, method, adata.n_obs)
            embeddings[method] = emb

            adata_tmp = adata.copy()
            sc.pp.neighbors(adata_tmp, n_pcs=min(50, X_pca.shape[1]), n_neighbors=N_NEIGHBORS)
            sc.tl.leiden(adata_tmp, resolution=0.5, random_state=RANDOM_STATE, flavor="igraph", n_iterations=2)
            leiden_labels = adata_tmp.obs["leiden"].astype(int).values

            ari = adjusted_rand_score(true_labels, leiden_labels)
            nmi = normalized_mutual_info_score(true_labels, leiden_labels)

            results["methods"][method] = {
                "runtime_s": round(runtime, 2),
                "ari": round(ari, 3),
                "nmi": round(nmi, 3),
            }
            print(f"    {runtime:6.2f}s | ARI: {ari:.3f} | NMI: {nmi:.3f}")

        # Save JSON
        with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
            json.dump(results, f, indent=2)

        # PLOTTING 
        n_types = results["n_types"]
        palette = "tab20" if n_types <= 20 else "hsv"
        show_legend = n_types <= 20

        # Dynamic figure size: wider if legend
        base_width = 15
        fig_width = base_width + 4 if show_legend else base_width
        fig_height = 6

        fig = plt.figure(figsize=(fig_width, fig_height))
        if show_legend:
            gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.35], wspace=0.35, hspace=0.3)
            axes = [fig.add_subplot(gs[i]) for i in range(3)]
            legend_ax = fig.add_subplot(gs[3])
            legend_ax.axis('off')
        else:
            gs = fig.add_gridspec(1, 3, wspace=0.35)
            axes = [fig.add_subplot(gs[i]) for i in range(3)]

        for idx, (m, emb) in enumerate(embeddings.items()):
            ax = axes[idx]
            df = pd.DataFrame({
                "x": emb[:, 0],
                "y": emb[:, 1],
                "cell_type": adata.obs["cell_type"].values
            })

            point_size = 8 if adata.n_obs > 10000 else 20
            sns.scatterplot(
                data=df, x="x", y="y",
                hue="cell_type" if show_legend else None,
                ax=ax,
                s=point_size,
                alpha=0.7,
                palette=palette,
                edgecolor="none",
                legend=show_legend
            )

            rt = results["methods"][m]["runtime_s"]
            ax.set_title(f"{m.upper()}\n({rt:.1f}s)", fontsize=12, pad=18)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis='both', which='major', labelsize=10)
            if show_legend:
                ax.get_legend().remove()

        # Unified legend on the right
        if show_legend:
            handles, labels = axes[0].get_legend_handles_labels()
            legend_ax.legend(
                handles, labels,
                loc='center left',
                title="cell_type",
                fontsize=9,
                title_fontsize=10,
                frameon=False,
                ncol=1,
                borderaxespad=0
            )

        # Main title with padding
        fig.suptitle(
            f"{ds_name.upper()} (n={adata.n_obs:,}, {n_types} types)",
            fontsize=14,
            y=0.98,
            weight='bold'
        )

        # Final layout
        plt.subplots_adjust(
            top=0.88,
            bottom=0.12,
            left=0.05,
            right=0.84 if show_legend else 0.95,
            hspace=0.3,
            wspace=0.35
        )

        plt.savefig(os.path.join(OUT_DIR, "embeddings.png"), dpi=300, bbox_inches="tight")
        plt.close()

        # APPEND ONLY ON SUCCESS 
        all_results.append(results)

    except Exception as e:
        print(f"  ERROR in {ds_name}: {e}")
        continue


# SCALING SUMMARY 
print(f"\n{'='*60}")
print("SCALING SUMMARY")
print(f"{'='*60}")

if not all_results:
    print("No results to summarize. Check errors above.")
else:
    summary = []
    for r in all_results:
        for m, v in r["methods"].items():
            summary.append({
                "Dataset": r["dataset"].upper(),
                "Method": m.upper(),
                "n_cells": r["n_cells"],
                "Runtime (s)": v["runtime_s"],
                "ARI": v["ari"],
                "NMI": v["nmi"],
            })

    df = pd.DataFrame(summary)
    print(df.to_string(index=False))
    df.to_csv("scaling_summary.csv", index=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    sns.lineplot(data=df, x="n_cells", y="ARI", hue="Method", marker="o", ax=ax1)
    ax1.set_title("ARI vs Dataset Size")
    ax1.set_xscale("log")

    sns.lineplot(data=df, x="n_cells", y="Runtime (s)", hue="Method", marker="o", ax=ax2)
    ax2.set_title("Runtime vs Dataset Size (log-log)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    plt.tight_layout()
    plt.savefig("scaling_plot.png", dpi=300, bbox_inches="tight")
    plt.close()

    print("\nSUCCESS! All outputs saved.")
    print("→ results_pbmc3k/ | results_pbmc68k/")
    print("→ scaling_summary.csv | scaling_plot.png")
