"""
Phase 3, Step 2: Cluster RBPs by ESM-2 embedding and propose phage cocktails.

WHY THIS STEP EXISTS
--------------------
Phage A and Phage B will infect overlapping strains of bacteria if their
Receptor Binding Proteins (RBPs) are similar — because they will both
recognize the same bacterial surface molecules. If their RBPs are
dissimilar, they recognize different bacterial features, so a cocktail
of both extends host range.

This step takes the ESM-2 embeddings from Step 1 and computes:

1. PAIRWISE PHAGE DISTANCE — for each pair of phages, the average
   minimum cosine distance from each RBP in phage A to the closest RBP
   in phage B. Low distance = redundant; high distance = complementary.

2. RBP CLUSTERING — group all RBPs across phages into clusters by
   embedding similarity. Each cluster represents a putative shared
   recognition mechanism.

3. COCKTAIL RECOMMENDATION — a greedy algorithm: start with the
   top-ranked phage from Phase 1, then iteratively add the phage that
   contributes the MOST NEW RBP CLUSTERS (i.e. the highest-coverage
   addition). Output 2-, 3-, and 4-phage cocktails.

OUTPUT INTERPRETATION
---------------------
A 2-phage cocktail of {A511, P100} that shares all RBP clusters is
worse than {A511, P70} where P70 brings novel RBPs the cocktail did
not previously cover — even though A511+P100 have higher individual
scores. This is the key trade-off for PHIND cartridge design.

Usage
-----
    python src/host_range/02_analyze_host_range.py

Inputs
------
    results/host_range/rbp_embeddings.npy   (from Step 1)
    results/host_range/rbp_metadata.csv     (from Step 1)
    results/candidate_phages_ranked.csv

Outputs
-------
    results/host_range/phage_distance_matrix.csv   — pairwise distance table
    results/host_range/rbp_clusters.csv            — per-RBP cluster assignments
    results/host_range/cocktail_recommendations.csv — proposed cocktails
    results/host_range/rbp_umap.png                — 2-D RBP map
    results/host_range/phage_distance_heatmap.png  — heatmap of phage pairs
    results/08_host_range_log.txt                  — summary
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_distances

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EMB_PATH = PROJECT_ROOT / "results" / "host_range" / "rbp_embeddings.npy"
META_PATH = PROJECT_ROOT / "results" / "host_range" / "rbp_metadata.csv"
RANKED_CSV = PROJECT_ROOT / "results" / "candidate_phages_ranked.csv"

OUT_DIR = PROJECT_ROOT / "results" / "host_range"
DIST_CSV = OUT_DIR / "phage_distance_matrix.csv"
CLUSTERS_CSV = OUT_DIR / "rbp_clusters.csv"
COCKTAIL_CSV = OUT_DIR / "cocktail_recommendations.csv"
UMAP_PNG = OUT_DIR / "rbp_umap.png"
HEATMAP_PNG = OUT_DIR / "phage_distance_heatmap.png"
LOG_OUT = PROJECT_ROOT / "results" / "08_host_range_log.txt"

# Brand
TEAL = "#5BA8A0"
TEAL_DARK = "#2D6E68"
NAVY = "#1A2F4D"
ACCENT = "#E6A85C"

# Number of RBP clusters to extract — chosen so each cluster has ~3-5 members
# on average for our 38-protein dataset. Tunable.
N_CLUSTERS = 8


# ---------------------------------------------------------------------------
# Core analyses
# ---------------------------------------------------------------------------

def compute_phage_distance(
    embeddings: np.ndarray, meta: pd.DataFrame,
) -> pd.DataFrame:
    """
    For each pair (A, B) of phages, average minimum cosine distance from
    each RBP in A to the closest RBP in B (and vice versa, symmetrized).
    """
    accessions = sorted(meta["accession"].unique())
    n = len(accessions)
    dist = np.zeros((n, n))

    for i, a in enumerate(accessions):
        emb_a = embeddings[meta["accession"] == a]
        for j, b in enumerate(accessions):
            if i == j:
                continue
            emb_b = embeddings[meta["accession"] == b]
            # Cosine distance between every (a, b) pair
            d = cosine_distances(emb_a, emb_b)
            # For each RBP in a, minimum distance to any RBP in b
            min_a_to_b = d.min(axis=1).mean()
            dist[i, j] = min_a_to_b

    # Symmetrize
    dist = (dist + dist.T) / 2
    return pd.DataFrame(dist, index=accessions, columns=accessions)


def cluster_rbps(embeddings: np.ndarray, n_clusters: int = N_CLUSTERS) -> np.ndarray:
    """Hierarchical agglomerative clustering of RBP embeddings."""
    clf = AgglomerativeClustering(n_clusters=n_clusters, metric="cosine", linkage="average")
    labels = clf.fit_predict(embeddings)
    return labels


def cocktail_greedy(
    meta: pd.DataFrame, ranked_order: list[str], max_size: int = 4,
) -> list[dict]:
    """
    Greedy cocktail expansion:
      - Seed cocktail with the top-1 phage from Phase 1 ranking
      - At each step add the phage that contributes the most NEW
        RBP clusters
    """
    cocktail = [ranked_order[0]]
    covered_clusters = set(
        meta[meta["accession"] == ranked_order[0]]["rbp_cluster"].tolist()
    )

    history = [
        {
            "size": 1,
            "phages": cocktail.copy(),
            "covered_clusters": len(covered_clusters),
            "added_phage": ranked_order[0],
            "new_clusters_added": len(covered_clusters),
        }
    ]

    for _ in range(max_size - 1):
        best_phage = None
        best_new = -1
        for cand in ranked_order:
            if cand in cocktail:
                continue
            cand_clusters = set(meta[meta["accession"] == cand]["rbp_cluster"].tolist())
            new = len(cand_clusters - covered_clusters)
            if new > best_new:
                best_new = new
                best_phage = cand
        if best_phage is None:
            break
        cocktail.append(best_phage)
        added_clusters = set(meta[meta["accession"] == best_phage]["rbp_cluster"].tolist())
        covered_clusters |= added_clusters
        history.append(
            {
                "size": len(cocktail),
                "phages": cocktail.copy(),
                "covered_clusters": len(covered_clusters),
                "added_phage": best_phage,
                "new_clusters_added": best_new,
            }
        )

    return history


# ---------------------------------------------------------------------------
# Visualizations
# ---------------------------------------------------------------------------

def plot_rbp_pca(embeddings: np.ndarray, meta: pd.DataFrame, out: Path) -> None:
    """2-D PCA scatter of RBP embeddings, colored by phage."""
    pca = PCA(n_components=2)
    proj = pca.fit_transform(embeddings)

    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.patch.set_facecolor("white")

    colors = [TEAL, NAVY, ACCENT, TEAL_DARK, "#A0D8D2"]
    accessions = sorted(meta["accession"].unique())
    for color, acc in zip(colors, accessions):
        mask = meta["accession"].values == acc
        n_rbps = mask.sum()
        # Get a friendly label
        label_full = meta[mask]["phage_name"].iloc[0]
        label = label_full.replace("Listeria phage ", "").split(",")[0]
        ax.scatter(
            proj[mask, 0], proj[mask, 1],
            s=120, alpha=0.75, color=color,
            label=f"{label}  (n={n_rbps})",
            edgecolors="white", linewidth=1.5,
        )

    ax.set_xlabel(f"PC1  ({pca.explained_variance_ratio_[0]*100:.1f}%)",
                  fontsize=11, color=NAVY)
    ax.set_ylabel(f"PC2  ({pca.explained_variance_ratio_[1]*100:.1f}%)",
                  fontsize=11, color=NAVY)
    ax.set_title(
        "ESM-2 embeddings of RBPs across top 5 candidates\n"
        "(closely-placed points = likely overlapping host range)",
        fontsize=13, color=NAVY, fontweight="bold", pad=12,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9, frameon=True)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close()


def plot_distance_heatmap(dist_df: pd.DataFrame, name_lookup: dict, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 6.5), dpi=150)
    fig.patch.set_facecolor("white")

    short_labels = [
        name_lookup.get(a, a).replace("Listeria phage ", "").split(",")[0][:18]
        for a in dist_df.index
    ]

    im = ax.imshow(dist_df.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(dist_df)))
    ax.set_yticks(range(len(dist_df)))
    ax.set_xticklabels(short_labels, rotation=35, ha="right", color=NAVY)
    ax.set_yticklabels(short_labels, color=NAVY)

    # Cell annotations
    for i in range(len(dist_df)):
        for j in range(len(dist_df)):
            color = "white" if dist_df.values[i, j] < dist_df.values.mean() else "black"
            ax.text(j, i, f"{dist_df.values[i, j]:.3f}",
                    ha="center", va="center", color=color, fontsize=9)

    plt.colorbar(im, ax=ax, label="Mean cosine distance")
    ax.set_title(
        "Pairwise RBP-repertoire distance between top 5 candidates\n"
        "(higher = more complementary for cocktail design)",
        fontsize=12, color=NAVY, fontweight="bold", pad=12,
    )
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    embeddings = np.load(EMB_PATH)
    meta = pd.read_csv(META_PATH)
    ranked = pd.read_csv(RANKED_CSV)

    print(f"\n[1/4] Loaded {embeddings.shape[0]} RBP embeddings "
          f"of dim {embeddings.shape[1]}.")

    # Pairwise phage distance
    print(f"\n[2/4] Computing pairwise phage RBP-repertoire distance...")
    dist_df = compute_phage_distance(embeddings, meta)
    dist_df.to_csv(DIST_CSV)
    print(f"      → wrote {DIST_CSV}")

    # RBP clustering
    print(f"\n[3/4] Clustering {len(embeddings)} RBPs into {N_CLUSTERS} groups...")
    labels = cluster_rbps(embeddings, N_CLUSTERS)
    meta = meta.copy()
    meta["rbp_cluster"] = labels
    meta.to_csv(CLUSTERS_CSV, index=False)
    print(f"      → wrote {CLUSTERS_CSV}")

    # Visualizations
    print(f"\n[4/4] Generating visualizations...")
    plot_rbp_pca(embeddings, meta, UMAP_PNG)
    name_lookup = dict(zip(ranked["accession"], ranked["name"]))
    plot_distance_heatmap(dist_df, name_lookup, HEATMAP_PNG)
    print(f"      → wrote {UMAP_PNG.name}, {HEATMAP_PNG.name}")

    # Cocktail recommendation
    ranked_order = ranked["accession"].tolist()
    # Filter to phages we actually embedded
    ranked_order = [a for a in ranked_order if a in meta["accession"].unique()]
    cocktail_history = cocktail_greedy(meta, ranked_order, max_size=len(ranked_order))

    cocktail_df = pd.DataFrame(cocktail_history)
    cocktail_df["phages"] = cocktail_df["phages"].apply(lambda lst: " + ".join(lst))
    cocktail_df.to_csv(COCKTAIL_CSV, index=False)

    # ---- Print + log summary ----
    log_lines = [
        "Phase 3 Step 2: Host-range analysis + cocktail recommendation",
        "=============================================================",
        f"RBP proteins analyzed:     {len(meta)}",
        f"Phages analyzed:           {meta['accession'].nunique()}",
        f"Embedding model:           ESM-2 (facebook/esm2_t12_35M_UR50D, 480-dim)",
        f"RBP clusters identified:   {N_CLUSTERS}",
        "",
        "Pairwise RBP-repertoire distance (cosine):",
    ]
    for line in dist_df.round(3).to_string().split("\n"):
        log_lines.append("  " + line)

    log_lines += ["", "RBPs per cluster:"]
    for cl in sorted(meta["rbp_cluster"].unique()):
        members = meta[meta["rbp_cluster"] == cl]
        log_lines.append(f"  Cluster {cl}: {len(members)} RBPs across "
                         f"{members['accession'].nunique()} phages")

    log_lines += ["", "Greedy cocktail recommendation:"]
    log_lines.append(f"{'Size':>4}  {'Phages':<60}  {'Clusters covered':>16}  {'New added':>10}")
    log_lines.append("-" * 100)
    for row in cocktail_history:
        phages_str = " + ".join([
            name_lookup.get(p, p).replace("Listeria phage ", "").split(",")[0]
            for p in row["phages"]
        ])
        log_lines.append(
            f"{row['size']:>4}  {phages_str[:60]:<60}  "
            f"{row['covered_clusters']:>16}  +{row['new_clusters_added']:>9}"
        )

    print("\n" + "\n".join(log_lines))
    with LOG_OUT.open("w") as fh:
        fh.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
