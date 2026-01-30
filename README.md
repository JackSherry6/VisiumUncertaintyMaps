# Overview:
This pipeline performs perturbation-based clustering on 10x visium spatial transcriptomics data to quantify cluster stability and uncertainty. It takes a standard space ranger output folder and builds a set of 
CSV/Parquet files and plots describing local/global co-clustering stability, perturbation metadata, and neighborhood-based uncertainty metrics.
- Currently works for 16um and 8um fresh-frozen bin sizes (still being optimized for 2um and FFPE data)

# About:

## Motivation:
One of my goals during the 2025/26 winter break was to get familiar with processing spatial transcriptomic data. As I tested spatial clustering parameters by running small perturbations and comparing the resulting maps side-by-side, the idea occured to me that the consistency of a spot’s neighborhood might be turned into a quantitative measure of confidence. In theory, if a spot repeatedly clusters with the same nearby spots across reasonable parameter variations, that suggests a stable signal, whereas if its assignments fluctuate, it may reflect boundary effects, mixed states, or technical noise.

## Problem:
Spot-based spatial transcriptomics trades resolution for scale, especially in-situ capturing (ISC) methods like 10x visium. While 10x has greatly shrunk the bin size all the way down from 16x16 to 2x2 um, the issue remains that bins can capture partial cells, overlapping cell processes, diffusion of transcripts, and segmentation or alignment noise, all of which contribute to some bins being an aggregate rather than a true single-cell readout, not to mention, the computational burdon for a 2x2 workflow can be extremely high. By standard workflows these mixed, ambiguous signals can be confidently misassigned to the wrong cell type or spatial domain. Deconvolution is commonly employed to fix this issue but it is not a perfect solution. As I found out during my testing, deconvolution depends greatly on the quality of the reference, can be sensitive to modeling assumptions, and doesn't do much to indicate which bins are robust versus ambiguous. This creates a practical issue, knowing where in the tissue the inferred labels are trustworthy, and where they should be treated cautiously—especially when downstream biological conclusions may hinge on small spatial boundaries or rare populations.

## Solution:
Attempt to create perturbations of various clustering parameters, track the cluster variation of each spot, and measure the consistency of each spot's neighborhood over perturbations using a kNN-type model. 

### Building the Perturbations (technical aspects):
- Perturbations are defined in the config.py file and built up-front in a single table (perturbations.parquet). Each row is one fully specified clustering run (unique run_id + parameter set). This makes the experiment reproducible and lets you submit the whole grid as an array job.
- The script takes a 1-based TASK_ID, loads exactly that row from perturbations.parquet, and runs it end-to-end. This helps with run time and makes it SGE/array-friendly.
- The specified parameters in config.py are intentionally “reasonable workflow variations,” but they can be changed based on situational preferences.
- Rather than hard thresholds that don’t scale across tissue/bin sizes, this pipeline computes per-run cutoffs via dataset quantiles of n_genes_by_counts and total_counts (gene_q, umi_q) instead of fixed cutoffs. This makes QC perturbations more comparable across datasets.
- I set min_cells = max(3, int(0.01 * n_obs)) so gene filtering scales with the number of spots and remains stable across different QC stringencies.
- Each run varies in n_hvgs which is capped by n_vars and attempts seurat_v3 with a fallback to cell_ranger. This had to be done to keep hvgs robust and prevent failure cases from biasing stability or completely breaking the clustering.
- When smooth > 0, the pipeline constructs a grid-based spatial adjacency, degree-normalizes it, and applies it iteratively "smooth times" so each spot is updated by a weighted average of its neighbors. This is done to prevent edge spots from being systematically biased because they have fewer neighbors than central spots.
- The pipeline switched scaling for memory-heavy perturbations (n_obs > 50,000) to sparse-friendly StandardScaler(with_mean=False) (otherwise sc.pp.scale) and clip values to [-10,10]. This allows my pipeline to run reasonably with larger 8um samples.
- n_pca varies per run but it caps at min(n_obs-1, n_vars-1) to avoid invalid PCA settings.
- The pipeline picks the PCA solver based on dataset size, trying randomized first for large datasets (otherwise arpack), and switching to the other solver if the first one fails.
- Random_seed is fixed to 42 for reproducibility.
- Each run writes either a success parquet or a failure parquet with the status, error, and parameters, and it uses .started and .tmp_failure files to flag incomplete runs. Runs that end up with too few spots or too few genes/HVGs after QC are flagged for error and stopped.
- Each perturbation writes a spot-level parquet with spot_id, cluster labels, cluster sizes, spatial coordinates, and the full parameter set to enable per-spot stability computation across runs. These parquets are merged into a single perurbations_expanded.parquet for downstream analysis.

### Calculating Stability (technical aspects: 

### Speed and Memory Optimizations:

## Example Figures:
-Soybean Data:
<table>
  <tr>
    <th align="center">Soybean 16 µm bins</th>
    <th align="center">Soybean 8 µm bins</th>
  </tr>
  <tr>
    <td align="center">
      <img src="example_plots/soybean_16um_spatial_scatterplot_ambiguity.png" width="350">
    </td>
    <td align="center">
      <img src="example_plots/soybean_8um_spatial_scatterplot_ambiguity.png" width="350">
    </td>
  </tr>
</table>

<table>
  <tr>
    <th align="center">Ovarian Cancer 16 µm bins</th>
    <th align="center">Ovarian Cancer 8 µm bins</th>
  </tr>
  <tr>
    <td align="center">
      <img src="example_plots/ovarian_cancer_16um_spatial_scatterplot_ambiguity.png" width="350">
    </td>
    <td align="center">
      <img src="example_plots/ovarian_cancer_8um_spatial_scatterplot_ambiguity.png" width="350">
    </td>
  </tr>
</table>
