# Diffusion pseudotime

The working analysis uses the existing MG-selected Seurat object, 20 PCs,
60 nearest neighbors, and 30 diffusion components. It adds no cell filtering.
The root is the medoid of the neighborhood maximizing the mean standardized MG
score minus the mean standardized neuronal-union score. The neuronal signature
is the union of the six curated retinal-neuron subtype sets. Progenitor and
neuronal neighborhoods are reference locations, not constrained endpoints.

`analysis.R` fits the trajectory and repeats root selection at 30 and 90
neighbors. These sensitivity runs change both graph size and the root. It saves
cell-level pseudotime, signatures, settings, diagnostics, and the diffusion fit,
plus a distance-colored UMAP and unbinned marker bars for the final cell quartile.
Pseudotime is divided by its maximum for the distance-colored plots; percentile
plots use the rank ordering instead. Marker bars show counts per 10,000 total RNA
counts, with one equal-width slot per cell.

`plot-marker-progression.R` reads that saved fit and produces 100 equal-cell bins
of mean normalized counts, including zeros, and detection fractions. Detection
means raw RNA count greater than zero. `plot-clusters.R` recreates the rank-colored
UMAP, configured cluster-colored UMAP, and cluster beeswarm without refitting.
The clustering is only a display annotation, not an input to DPT.

Run the scripts with `source()` from the project root, with `analysis.R` first.
The plotting scripts do not refit DPT. The exploratory dependencies are listed
in `DESCRIPTION` under `Suggests`.

All scripts default to `outputs/mg-minus-neuron/`. To run without touching saved
results, set `ESPI_PSEUDOTIME_OUTPUT_DIR` to a separate directory and run the
analysis followed by the plotting scripts. Existing files still require explicit
`ESPI_OVERWRITE=true`; publication outputs are not used.

Root-level `outputs/` preserves the earlier MG-minus-progenitor run and its
interactive plotting variants. Existing `marker_trends.csv` snapshots use the
original 20 log-expression bins; current fitting code no longer exports those
summaries. The two original main runs selected the same root and cell ordering.
Older snapshots are retained rather than rewritten to match later plot styling.
