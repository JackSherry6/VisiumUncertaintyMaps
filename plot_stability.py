#!/usr/bin/env python

import argparse
import os
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N_COLS = 6
TILE_INCH = 2.1
DPI = 250
INVERT_Y = True
CMAP = "viridis"
NORMALIZE_WITHIN_CLUSTER = True
BG_S, BG_ALPHA = 0.6, 0.05
CLUSTER_S, CLUSTER_ALPHA = 0.5, 0.95

def parse_args():
    p = argparse.ArgumentParser(description="Plot spatial stability/ambiguity metrics.")
    p.add_argument("--stab_csv", default="spot_jaccard_stability.csv", type=str)
    p.add_argument("--pert_parquet", default="perturbations_expanded.parquet", type=str)
    p.add_argument("--outdir", default="figures", type=str)
    p.add_argument("--bins", default=50, type=int)
    p.add_argument("--point_size", default=0.3, type=float)
    return p.parse_args()

def metric_specs(stab_cols):
    specs = []

    if "ambiguity" in stab_cols:
        specs.append(
            ("ambiguity", "Ambiguity (higher = worse)",
             "Spot ambiguity distribution", "Spatial spot ambiguity")
        )
    if "local_stability" in stab_cols:
        specs.append(
            ("local_stability", "Local co-clustering stability (higher = better)",
             "Local stability distribution", "Spatial local stability")
        )
    if "global_stability" in stab_cols:
        specs.append(
            ("global_stability", "Global co-clustering stability (higher = better)",
             "Global stability distribution", "Spatial global stability")
        )
    if "expected_jaccard_stability" in stab_cols:
        specs.append(
            ("expected_jaccard_stability", "Co-clustering stability (higher = better)",
             "Spot stability distribution", "Spatial spot stability")
        )

    return specs

def clusters_csv_from_config() -> Path:
    try:
        import config
    except Exception as e:
        raise ImportError(f"Could not import config.py ({e}).")

    cfg = getattr(config, "CONFIG", None)
    if not isinstance(cfg, dict):
        raise ValueError("config.py must define a dict named CONFIG.")

    try:
        data_path = cfg["paths"]["data_path"]
    except Exception:
        raise KeyError("CONFIG must contain CONFIG['paths']['data_path'].")

    clusters_csv = (
        Path(data_path)
        / "analysis"
        / "clustering"
        / "gene_expression_graphclust"
        / "clusters.csv"
    )
    return clusters_csv

def load_clusters_raw(clusters_csv: Path) -> pd.DataFrame:
    cl = pd.read_csv(clusters_csv, compression="infer")

    if "Barcode" in cl.columns and "Cluster" in cl.columns:
        cl = cl.rename(columns={"Barcode": "spot_id", "Cluster": "cluster"})
    elif "spot_id" in cl.columns and "cluster" in cl.columns:
        pass
    else:
        raise ValueError(
            f"{clusters_csv} must contain either (Barcode, Cluster) or (spot_id, cluster).\n"
            f"Columns found: {list(cl.columns)}"
        )

    cl = cl[["spot_id", "cluster"]].copy()
    cl["spot_id"] = cl["spot_id"].astype(str)

    cl["cluster"] = (
        cl["cluster"]
        .astype(str)
        .str.extract(r"(-?\d+)", expand=False)
    )
    cl["cluster"] = pd.to_numeric(cl["cluster"], errors="coerce")
    return cl


def harmonize_cluster_spot_ids(clusters: pd.DataFrame, coords_spot_ids: pd.Series) -> pd.DataFrame:
    coords_ids = set(coords_spot_ids.astype(str).tolist())
    s = clusters["spot_id"].astype(str)

    def score(series: pd.Series) -> int:
        return int(series.isin(coords_ids).sum())

    cand = []
    cand.append(("identity", s))
    cand.append(("add_-1", s.where(s.str.contains(r"-\d+$"), s + "-1")))
    cand.append(("strip_-digits", s.str.replace(r"-\d+$", "", regex=True)))
    cand.append(("strip_-digits_then_add_-1", s.str.replace(r"-\d+$", "", regex=True) + "-1"))

    scored = [(name, score(v)) for name, v in cand]
    best_name, best_score = max(scored, key=lambda x: x[1])

    clusters = clusters.copy()
    clusters["spot_id"] = dict(cand)[best_name]

    total = len(clusters)
    print(f"[clusters] spot_id harmonization: {best_name} (overlap {best_score}/{total})")
    return clusters

def make_grid(n_tiles: int, n_cols: int):
    n_rows = math.ceil(n_tiles / n_cols)
    fig_w = n_cols * TILE_INCH
    fig_h = n_rows * TILE_INCH
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), constrained_layout=True)
    axes = np.array(axes).reshape(-1)
    return fig, axes


def style_tile(ax):
    ax.set_aspect("equal")
    if INVERT_Y:
        ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    stab = pd.read_csv(args.stab_csv)
    if "spot_id" not in stab.columns:
        raise ValueError(f"{args.stab_csv} must contain a 'spot_id' column.")
    stab["spot_id"] = stab["spot_id"].astype(str)

    specs = metric_specs(stab.columns)
    if not specs:
        raise KeyError(
            "Could not find any expected stability columns in stability CSV. "
            "Looked for: ambiguity, local_stability, global_stability, expected_jaccard_stability."
        )

    print("Will plot metrics:", ", ".join([c for c, *_ in specs]))

    coords = pd.read_parquet(args.pert_parquet)[["spot_id", "x", "y"]].drop_duplicates("spot_id").copy()
    coords["spot_id"] = coords["spot_id"].astype(str)
    
    plot_df = stab.merge(coords, on="spot_id", how="inner")
    if plot_df.empty:
        raise RuntimeError("No overlapping spot_id between stability and coordinates")
        
    for metric_col, metric_label, title_hist, title_spatial in specs:
        series = plot_df[metric_col].dropna()
        if series.empty:
            print(f"Skipping '{metric_col}' (all values are NaN after merge).")
            continue

        print(f"Plotting '{metric_col}' (n={series.shape[0]}) ...")

        plt.figure(figsize=(4, 3))
        plt.hist(series, bins=args.bins)
        plt.xlabel(metric_label)
        plt.ylabel("Number of spots")
        plt.title(title_hist)
        plt.tight_layout()
        plt.savefig(Path(args.outdir) / f"stability_histogram_{metric_col}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(6, 6))
        scat = plt.scatter(plot_df["x"], plot_df["y"], c=plot_df[metric_col], s=args.point_size,)
        plt.colorbar(scat, label=metric_label)
        plt.gca().invert_yaxis()
        plt.axis("equal")
        plt.axis("off")
        plt.title(title_spatial)
        plt.tight_layout()
        plt.savefig(Path(args.outdir) / f"spatial_scatterplot_{metric_col}.png", dpi=300)
        plt.close()

    print(f"Done. Plotted {len(plot_df)} spots (after merge).")
    
    if "ambiguity" not in plot_df.columns:
        raise ValueError(
            "Cannot create ambiguity gallery because 'ambiguity' column is missing from stability CSV."
        )

    clusters_path = clusters_csv_from_config()
    if not clusters_path.exists():
        raise FileNotFoundError(
            "clusters.csv not found at the config-derived path:\n"
            f"  {clusters_path}\n"
            "Check CONFIG['paths']['data_path'] in config.py and confirm clustering output exists there."
        )

    print(f"[gallery] Using clusters file (from config.py): {clusters_path}")
    clusters = load_clusters_raw(clusters_path)
    clusters = harmonize_cluster_spot_ids(clusters, coords["spot_id"])

    spot_df = coords.merge(clusters, on="spot_id", how="left").merge(
        plot_df[["spot_id", "ambiguity"]], on="spot_id", how="left"
    )

    n_with_cluster = int(spot_df["cluster"].notna().sum())
    if n_with_cluster == 0:
        coords_head = coords["spot_id"].head(5).tolist()
        cl_head = clusters["spot_id"].head(5).tolist()
        raise ValueError(
            "[gallery] clusters.csv loaded, but 0 spots matched after merge.\n"
            f"  coords spot_id examples:   {coords_head}\n"
            f"  clusters spot_id examples: {cl_head}\n"
            "This usually means clusters.csv is from a different run/dataset, or barcodes still differ."
        )

    cluster_ids = spot_df["cluster"].dropna().astype(int).sort_values().unique().tolist()
    if not cluster_ids:
        raise ValueError("[gallery] No clusters found (or all clusters parsed as NaN).")

    print(f"[gallery] Found {len(cluster_ids)} clusters (matched clusters for {n_with_cluster} spots)")

    bg = spot_df.dropna(subset=["x", "y"])
    bgx = bg["x"].to_numpy()
    bgy = bg["y"].to_numpy()

    fig_a, axes_a = make_grid(len(cluster_ids), N_COLS)
    fig_a.suptitle("Cluster-specific ambiguity", fontsize=14)

    last_scatter = None

    for i, cl in enumerate(cluster_ids):
        ax = axes_a[i]

        cl_spots = spot_df[spot_df["cluster"] == cl].dropna(subset=["x", "y"])
        ax.scatter(bgx, bgy, s=BG_S, alpha=BG_ALPHA, color="gray")

        cl_amb = cl_spots.dropna(subset=["ambiguity"]).copy()
        if len(cl_amb) == 0:
            ax.set_title(f"{cl} (no amb)", fontsize=9)
            style_tile(ax)
            continue

        if NORMALIZE_WITHIN_CLUSTER:
            amin = cl_amb["ambiguity"].min()
            amax = cl_amb["ambiguity"].max()
            if amax > amin:
                cvals = (cl_amb["ambiguity"] - amin) / (amax - amin)
            else:
                cvals = np.full(len(cl_amb), 0.5)
            vmin, vmax = 0.0, 1.0
        else:
            cvals = cl_amb["ambiguity"].to_numpy()
            vmin = float(spot_df["ambiguity"].min())
            vmax = float(spot_df["ambiguity"].max())

        last_scatter = ax.scatter(
            cl_amb["x"], cl_amb["y"],
            s=CLUSTER_S, alpha=CLUSTER_ALPHA,
            c=cvals, cmap=CMAP, vmin=vmin, vmax=vmax
        )

        ax.set_title(f"{cl} (n={len(cl_amb)})", fontsize=9)
        style_tile(ax)

    for j in range(len(cluster_ids), len(axes_a)):
        axes_a[j].axis("off")

    if last_scatter is not None:
        cbar = fig_a.colorbar(last_scatter, ax=axes_a[:len(cluster_ids)], shrink=0.85, pad=0.01)
        cbar.set_label("ambiguity (cluster-normalized)" if NORMALIZE_WITHIN_CLUSTER else "ambiguity (global scale)")

    out_path = Path(args.outdir) / "ambiguity_gallery.png"
    fig_a.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig_a)
    print(f"[gallery] Saved {out_path}")

    print("All done.")


if __name__ == "__main__":
    main()
