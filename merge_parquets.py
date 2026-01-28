import pandas as pd
from pathlib import Path

pert = pd.read_parquet("perturbations.parquet")

dfs = []
for p in Path("perturbation_outputs").glob("*.parquet"):
    df = pd.read_parquet(p)
    assert df["run_id"].nunique() == 1, f"Multiple run_ids in {p}"
    dfs.append(df)

if not dfs:
    raise RuntimeError("No perturbation output files found.")

merged = pd.concat(dfs, ignore_index=True)

missing = set(pert["run_id"]) - set(merged["run_id"])
if missing:
    raise RuntimeError(f"Missing perturbations: {missing}")

merged.to_parquet("perturbations_expanded.parquet", index=False)

out_ids = merged["run_id"].unique()
print("All perturbations completed:", set(out_ids) == set(pert["run_id"]))
print("Merged shape:", merged.shape)
