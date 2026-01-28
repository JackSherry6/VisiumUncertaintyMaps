#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os

os.makedirs("figures", exist_ok=True)
stab = pd.read_csv("spot_jaccard_stability.csv")

metric_col = None
metric_label = None
plot_title_hist = None
plot_title_spatial = None

if "ambiguity" in stab.columns:
    metric_col = "ambiguity"
    metric_label = "Ambiguity (1 - stability)"
    plot_title_hist = "Spot ambiguity distribution"
    plot_title_spatial = "Spatial spot ambiguity"
elif "local_stability" in stab.columns:
    metric_col = "local_stability"
    metric_label = "Local co-clustering stability"
    plot_title_hist = "Local stability distribution"
    plot_title_spatial = "Spatial local stability"
elif "global_stability" in stab.columns:
    metric_col = "global_stability"
    metric_label = "Global co-clustering stability"
    plot_title_hist = "Global stability distribution"
    plot_title_spatial = "Spatial global stability"
elif "expected_jaccard_stability" in stab.columns:
    # Backwards-compatible with old outputs
    metric_col = "expected_jaccard_stability"
    metric_label = "Co-clustering stability"
    plot_title_hist = "Spot stability distribution"
    plot_title_spatial = "Spatial spot stability"
else:
    raise KeyError(
        "Could not find any of the expected stability columns in "
        "spot_jaccard_stability.csv. "
        "Looked for: ambiguity, local_stability, global_stability, "
        "expected_jaccard_stability."
    )

print(f"Using column '{metric_col}' for plotting.")

plt.figure(figsize=(4, 3))
plt.hist(stab[metric_col].dropna(), bins=50)
plt.xlabel(metric_label)
plt.ylabel("Number of spots")
plt.title(plot_title_hist)
plt.tight_layout()
plt.savefig("figures/stability_histogram.png", dpi=300)
plt.close()

coords = (
    pd.read_parquet("perturbations_expanded.parquet")
      [["spot_id", "x", "y"]]
      .drop_duplicates("spot_id")
)

plot_df = stab.merge(coords, on="spot_id", how="inner")

if plot_df.empty:
    raise RuntimeError("No overlapping spot_id between stability and coordinates")

plt.figure(figsize=(6, 6))
scat = plt.scatter(
    plot_df["x"],
    plot_df["y"],
    c=plot_df[metric_col],
    s=0.3,
)

cbar = plt.colorbar(scat, label=metric_label)
plt.gca().invert_yaxis()
plt.axis("equal")
plt.axis("off")
plt.title(plot_title_spatial)
plt.tight_layout()
plt.savefig("figures/spatial_scatterplot.png", dpi=300)
plt.close()

print(f"Plotted {len(plot_df)} spots using '{metric_col}'")
