#!/usr/bin/env python

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import json
import sys
import traceback
import warnings
from pathlib import Path
from scipy.sparse import csr_matrix, issparse
from sklearn.preprocessing import StandardScaler
from config import CONFIG

# suppress warnings for cleaner logs (just this once)
warnings.filterwarnings("ignore")

TASK_ID = int(sys.argv[1])
PERT_FILE = Path("perturbations.parquet")
OUTDIR = Path("perturbation_outputs")
OUTDIR.mkdir(exist_ok=True)

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

MIN_SPOTS_REQUIRED = 20
LARGE_NOBS_THRESHOLD = 50000  # when to switch to sparse-friendly scaling

def save_output(run_id, data_dict):
    try:
        df = pd.DataFrame(data_dict)
        df.to_parquet(OUTDIR / f"{run_id}.parquet", index=False)
        return True
    except Exception as e:
        print(f"ERROR saving output for {run_id}: {e}", file=sys.stderr)
        return False


def save_failure(run_id, status, error_msg, params, config):
    output = {
        "run_id": [run_id],
        "status": [status],
        "error": [str(error_msg)[:500]],
        **{k: [v] for k, v in params.items()},
        **{k: [v] for k, v in config["sample"].items()},
    }
    return save_output(run_id, output)

def load_visium(cfg):
    path = Path(cfg["paths"]["data_path"])
    lib = cfg["sample"]["library_id"]

    ad = sc.read_10x_h5(path / "filtered_feature_bc_matrix.h5")
    ad.var_names_make_unique()

    if not issparse(ad.X):
        ad.X = csr_matrix(ad.X)
    ad.X = ad.X.astype(np.float32)

    tissue = (
        pd.read_parquet(path / "spatial/tissue_positions.parquet")
        .set_index("barcode")
    )
    ad.obs = ad.obs.join(tissue, how="left")
    ad = ad[ad.obs["in_tissue"] == 1].copy()

    ad.obsm["spatial"] = ad.obs[
        ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    ].to_numpy()

    ad.uns["spatial"] = {lib: {}}
    with open(path / "spatial/scalefactors_json.json") as f:
        ad.uns["spatial"][lib]["scalefactors"] = json.load(f)

    return ad

#For smaller datasets: use scanpy's sc.pp.scale
#For large datasets: StandardScaler(with_mean=False) on sparse matrix
def scale_expression(adata):
    try:
        if adata.n_obs > LARGE_NOBS_THRESHOLD:
            # Sparse-friendly global scaling
            scaler = StandardScaler(with_mean=False, with_std=True)
            X_scaled = scaler.fit_transform(adata.X)

            if issparse(X_scaled):
                X_scaled.data = np.clip(X_scaled.data, -10, 10)
                adata.X = X_scaled.tocsr().astype(np.float32)
            else:
                X_scaled = np.clip(X_scaled, -10, 10).astype(np.float32)
                adata.X = csr_matrix(X_scaled)
        else:
            # Standard dense scaling for moderate-size datasets
            sc.pp.scale(adata, max_value=10)
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X.astype(np.float32))
            else:
                adata.X = adata.X.astype(np.float32)
    except Exception as e:
        raise RuntimeError(f"Scaling failed: {e}")

# Main processing:
def run_perturbation(params, run_id):
    heartbeat = OUTDIR / f"{run_id}.started"
    tmp_failure = OUTDIR / f"{run_id}.tmp_failure"
    heartbeat.touch()
    tmp_failure.touch()

    try:
        # Load data
        adata = load_visium(CONFIG)
        adata.uns["was_subsampled"] = False  # just for if we add subsampling later

        sc.pp.calculate_qc_metrics(adata, inplace=True)

        gene_cut = np.quantile(adata.obs.n_genes_by_counts, params["gene_q"])
        umi_cut = np.quantile(adata.obs.total_counts, params["umi_q"])
        mask = (
            (adata.obs.n_genes_by_counts >= gene_cut)
            & (adata.obs.total_counts >= umi_cut)
        )
        adata = adata[mask].copy()

        if adata.n_obs < MIN_SPOTS_REQUIRED:
            save_failure(
                run_id,
                "failed_qc_spots",
                f"Only {adata.n_obs} spots after QC",
                params,
                CONFIG,
            )
            return

        min_cells = max(3, int(0.01 * adata.n_obs))
        sc.pp.filter_genes(adata, min_cells=min_cells)

        if adata.n_vars < 50:
            save_failure(
                run_id,
                "failed_qc_genes",
                f"Only {adata.n_vars} genes after filtering",
                params,
                CONFIG,
            )
            return

        if params["norm"] == "log1p":
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        else:
            save_failure(
                run_id,
                "failed_invalid_norm",
                f"Unknown norm mode: {params['norm']}",
                params,
                CONFIG,
            )
            return

        if not issparse(adata.X):
            adata.X = csr_matrix(adata.X)
        adata.X = adata.X.astype(np.float32)

        n_top_genes = min(params.get("n_hvgs", 2000), adata.n_vars)

        hvg_success = False
        for flavor in ["seurat_v3", "cell_ranger"]:
            try:
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_top_genes,
                    flavor=flavor,
                    subset=True,
                )
                hvg_success = True
                break
            except Exception as e:
                if flavor == "cell_ranger":
                    save_failure(
                        run_id,
                        "failed_hvg_selection",
                        f"HVG selection failed: {str(e)[:200]}",
                        params,
                        CONFIG,
                    )
                    return

        if not hvg_success or adata.n_vars < 10:
            save_failure(
                run_id,
                "failed_hvg_count",
                f"Only {adata.n_vars} HVGs",
                params,
                CONFIG,
            )
            return

        if params.get("smooth", 0) > 0:
            try:
                sq.gr.spatial_neighbors(adata, coord_type="grid")
                adj = adata.obsp["spatial_connectivities"].astype(np.float32).tocsr()

                # Degree-normalized adjacency
                deg = np.asarray(adj.sum(axis=1)).flatten()
                deg[deg == 0] = 1.0
                adj = adj.multiply(1.0 / deg[:, None])

                for _ in range(params["smooth"]):
                    adata.X = (adj @ adata.X).tocsr().astype(np.float32)
            except Exception as e:
                save_failure(
                    run_id,
                    "failed_smoothing",
                    f"Smoothing failed: {str(e)[:200]}",
                    params,
                    CONFIG,
                )
                return

        try:
            scale_expression(adata)
        except RuntimeError as e:
            save_failure(run_id, "failed_scaling", str(e), params, CONFIG)
            return

        n_comps = min(
            params.get("n_pca", 30),
            adata.n_obs - 1,
            adata.n_vars - 1,
        )
        if n_comps < 2:
            save_failure(
                run_id,
                "failed_pca_setup",
                f"Invalid n_comps={n_comps}",
                params,
                CONFIG,
            )
            return

        solvers = (
            ["randomized", "arpack"]
            if adata.n_obs > LARGE_NOBS_THRESHOLD
            else ["arpack", "randomized"]
        )

        pca_success = False
        for solver in solvers:
            try:
                # Ensure sparse format for PCA
                if not issparse(adata.X):
                    adata.X = csr_matrix(adata.X)
                sc.pp.pca(
                    adata,
                    n_comps=n_comps,
                    svd_solver=solver,
                    random_state=RANDOM_SEED,
                )
                pca_success = True
                break
            except Exception as e:
                last_err = e

        if not pca_success:
            save_failure(
                run_id,
                "failed_pca",
                f"PCA failed with solvers {solvers}: {str(last_err)[:200]}",
                params,
                CONFIG,
            )
            return

        try:
            sc.pp.neighbors(
                adata,
                n_neighbors=params["n_neighbors"],
                random_state=RANDOM_SEED,
            )
            sc.tl.leiden(
                adata,
                resolution=params["leiden_res"],
                random_state=RANDOM_SEED,
            )
        except Exception as e:
            save_failure(
                run_id,
                "failed_clustering",
                f"Clustering failed: {str(e)[:200]}",
                params,
                CONFIG,
            )
            return

        if "leiden" not in adata.obs:
            save_failure(
                run_id,
                "failed_clustering",
                "Leiden produced no labels",
                params,
                CONFIG,
            )
            return

        sizes = adata.obs["leiden"].value_counts()
        n = adata.n_obs

        output_data = {
            "run_id": [run_id] * n,
            "spot_id": adata.obs_names.tolist(),
            "cluster": adata.obs["leiden"].astype(str).tolist(),
            "cluster_size": adata.obs["leiden"].map(sizes).tolist(),
            "x": adata.obsm["spatial"][:, 0].tolist(),
            "y": adata.obsm["spatial"][:, 1].tolist(),
            "status": ["success"] * n,
            "error": [""] * n,
            **{k: [v] * n for k, v in params.items()},
            **{k: [v] * n for k, v in CONFIG["sample"].items()},
            "subsampled": [adata.uns.get("was_subsampled", False)] * n,
        }

        tmp_failure.unlink(missing_ok=True)
        if not save_output(run_id, output_data):
            print(
                f"CRITICAL: Could not save success output for {run_id}",
                file=sys.stderr,
            )
            sys.exit(1)

    except Exception as e:
        # Catch any unexpected error during the main pipeline
        error_msg = f"{type(e).__name__}: {e}"
        traceback.print_exc()
        save_failure(run_id, "failed_exception", error_msg, params, CONFIG)

    finally:
        heartbeat.unlink(missing_ok=True)


if __name__ == "__main__":
    run_id = None
    try:
        df = pd.read_parquet(PERT_FILE)
        if TASK_ID < 1 or TASK_ID > len(df):
            raise IndexError(f"TASK_ID {TASK_ID} out of range (1..{len(df)})")

        params = df.iloc[TASK_ID - 1].to_dict()
        run_id = params.pop("run_id")

        run_perturbation(params, run_id)

    except Exception as e:
        print(f"CRITICAL ERROR in task {TASK_ID}: {e}", file=sys.stderr)
        traceback.print_exc()

        emergency_id = run_id if run_id else f"emergency_task_{TASK_ID}"
        try:
            emergency_output = {
                "run_id": [emergency_id],
                "status": ["failed_critical"],
                "error": [f"{type(e).__name__}: {str(e)}"[:500]],
                "task_id": [TASK_ID],
                **{k: [v] for k, v in CONFIG["sample"].items()},
            }
            save_output(emergency_id, emergency_output)
        except Exception:
            # If it gets to this point its completely doomed
            pass

        sys.exit(1)
