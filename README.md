## Overview:
This pipeline performs perturbation-based clustering on 10x visium spatial transcriptomics data to quantify cluster stability and uncertainty. It takes a standard space ranger output folder and builds a set of 
CSV/Parquet files and plots describing local/global co-clustering stability, perturbation metadata, and neighborhood-based uncertainty metrics.
- Currently works for 16um and 8um fresh-frozen bin sizes (still being optimized for 2um and FFPE data)

## Problem:
Spatial cluster assignment and deconvolution is...

## Solution:
Use spot stability scores to make deconvolution more robust
