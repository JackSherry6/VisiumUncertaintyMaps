#!/usr/bin/env python3
import argparse
import random
from itertools import combinations
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute global + local Jaccard stability and ambiguity."
    )
    parser.add_argument(
        "input_parquet",
        type=str,
        help="Merged perturbations parquet (from merge_parquets.qsub).",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        default="spot_jaccard_stability.csv",
        help="Output CSV file (default: spot_jaccard_stability.csv)",
    )
    parser.add_argument(
        "--max_run_pairs",
        type=int,
        default=200,
        help="Maximum number of run pairs to sample for stability computation.",
    )
    parser.add_argument(
        "--k_neighbors",
        type=int,
        default=6,
        help="Number of spatial neighbors per spot for local stability.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for run pair sampling.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.7,
        help="Weight for local uncertainty (1 - local_stability) in ambiguity.",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=0.3,
        help="Weight for global uncertainty (1 - global_stability) in ambiguity.",
    )
    return parser.parse_args()

def build_neighbors(coords: np.ndarray, k: int) -> np.ndarray:
    n_spots = coords.shape[0]
    n_neighbors = min(k + 1, n_spots)
    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto")
    nn.fit(coords)
    distances, indices = nn.kneighbors(coords)
    if n_neighbors > 1:
        neighbors = indices[:, 1:]
    else:
        neighbors = np.empty((n_spots, 0), dtype=np.int32)
    return neighbors.astype(np.int32)


def compute_global_stability(run_ids, run_to_clusters, n_spots: int, run_pairs,):
    global_sum = np.zeros(n_spots, dtype=np.float64)
    global_count = np.zeros(n_spots, dtype=np.int32)

    for idx, (ri, rj) in enumerate(run_pairs, 1):
        r1 = run_ids[ri]
        r2 = run_ids[rj]
        print(f"[global] Run pair {idx}/{len(run_pairs)}: {r1} vs {r2}")

        clust1 = run_to_clusters[ri]
        clust2 = run_to_clusters[rj]

        for c1_spots in clust1.values():
            len1 = c1_spots.size
            if len1 == 0:
                continue
            for c2_spots in clust2.values():
                len2 = c2_spots.size
                if len2 == 0:
                    continue

                inter_spots = np.intersect1d(
                    c1_spots, c2_spots, assume_unique=True
                )
                inter_len = inter_spots.size
                if inter_len == 0:
                    continue

                union = len1 + len2 - inter_len
                j = inter_len / union if union > 0 else 0.0
                if j == 0.0:
                    continue

                global_sum[inter_spots] += j
                global_count[inter_spots] += 1

    global_stability = np.zeros(n_spots, dtype=np.float64)
    mask = global_count > 0
    global_stability[mask] = global_sum[mask] / global_count[mask]

    return global_stability, global_count


def compute_local_stability(run_ids, labels, neighbors: np.ndarray, run_pairs,):
    n_runs, n_spots = labels.shape
    k = neighbors.shape[1]

    local_sum = np.zeros(n_spots, dtype=np.float64)
    local_count = np.zeros(n_spots, dtype=np.int32)

    for idx, (ri, rj) in enumerate(run_pairs, 1):
        r1 = run_ids[ri]
        r2 = run_ids[rj]
        print(f"[local] Run pair {idx}/{len(run_pairs)}: {r1} vs {r2}")

        labels_r1 = labels[ri]
        labels_r2 = labels[rj]

        for s in range(n_spots):
            nbr_idx = neighbors[s]  # (k,)

            if nbr_idx.size == 0:
                continue

            c1 = labels_r1[s]
            c2 = labels_r2[s]
            if c1 < 0 or c2 < 0:
                continue

            nbr_lab1 = labels_r1[nbr_idx]
            nbr_lab2 = labels_r2[nbr_idx]

            mask1 = nbr_lab1 == c1
            mask2 = nbr_lab2 == c2

            if not (mask1.any() or mask2.any()):
                continue

            inter = np.logical_and(mask1, mask2).sum()
            union = np.logical_or(mask1, mask2).sum()
            if union == 0:
                continue

            j = inter / union
            if j == 0.0:
                continue

            local_sum[s] += j
            local_count[s] += 1

    local_stability = np.zeros(n_spots, dtype=np.float64)
    mask = local_count > 0
    local_stability[mask] = local_sum[mask] / local_count[mask]

    return local_stability, local_count


def main():
    args = parse_args()

    print(f"Loading {args.input_parquet} ...")
    df = pd.read_parquet(
        args.input_parquet,
        columns=["run_id", "spot_id", "cluster", "x", "y"],
    )

    required_cols = {"run_id", "spot_id", "cluster", "x", "y"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Input parquet is missing required columns: {missing}. "
            "Make sure merge_parquets.py writes x, y, run_id, spot_id, cluster."
        )

    spot_codes, spot_ids = pd.factorize(df["spot_id"])
    df["spot_idx"] = spot_codes.astype(np.int32)
    n_spots = len(spot_ids)
    print(f"Total spots: {n_spots}")

    run_codes, run_ids = pd.factorize(df["run_id"])
    df["run_idx"] = run_codes.astype(np.int32)
    n_runs = len(run_ids)
    print(f"Total runs: {n_runs}")

    coords = np.zeros((n_spots, 2), dtype=np.float64)
    coords[df["spot_idx"].values, 0] = df["x"].values
    coords[df["spot_idx"].values, 1] = df["y"].values

    print(f"Building k-NN neighbors (k={args.k_neighbors}) ...")
    neighbors = build_neighbors(coords, k=args.k_neighbors)
    print("Neighbors built.")

    labels = np.full((n_runs, n_spots), -1, dtype=np.int32)
    run_to_clusters = [dict() for _ in range(n_runs)]

    print("Building label matrix and cluster→spots mapping ...")
    for r_idx, run in enumerate(run_ids):
        g = df[df["run_idx"] == r_idx]

        clust_codes, clust_labels = pd.factorize(g["cluster"])
        spot_idx = g["spot_idx"].values

        clust_codes = clust_codes.astype(np.int32)
        labels[r_idx, spot_idx] = clust_codes

        clusters = {}
        for code in np.unique(clust_codes):
            mask = clust_codes == code
            spots = spot_idx[mask]
            clusters[int(code)] = np.sort(spots.astype(np.uint32))
        run_to_clusters[r_idx] = clusters

    del df

    all_pairs = list(combinations(range(n_runs), 2))
    if len(all_pairs) == 0:
        raise ValueError("Need at least 2 runs to compute stability.")

    if len(all_pairs) > args.max_run_pairs:
        rng = random.Random(args.seed)
        run_pairs = rng.sample(all_pairs, args.max_run_pairs)
    else:
        run_pairs = all_pairs

    print(f"Using {len(run_pairs)} run pairs (out of {len(all_pairs)} total).")

    print("Computing GLOBAL cluster-based stability ...")
    global_stability, global_count = compute_global_stability(
        run_ids, run_to_clusters, n_spots, run_pairs
    )

    print("Computing LOCAL neighbor-based stability ...")
    local_stability, local_count = compute_local_stability(
        run_ids, labels, neighbors, run_pairs
    )

    alpha = args.alpha
    beta = args.beta
    ambiguity = alpha * (1.0 - local_stability) + beta * (1.0 - global_stability)

    out = pd.DataFrame(
        {
            "spot_id": spot_ids,
            "global_stability": global_stability,
            "local_stability": local_stability,
            "ambiguity": ambiguity,
            "n_pairs_global": global_count,
            "n_pairs_local": local_count,
        }
    )

    out.to_csv(args.output_csv, index=False)
    print(f"Saved stability + ambiguity to {args.output_csv}")


if __name__ == "__main__":
    main()
