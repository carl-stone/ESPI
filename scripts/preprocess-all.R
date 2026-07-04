#!/usr/bin/env Rscript

# Run every current preprocessing branch.
#
# Usage:
#   Rscript scripts/preprocess-all.R
#
# Arguments:
#   None.
#
# Branches run:
#   log1p without cell-cycle-HVG filtering
#   log1p with cell-cycle-HVG filtering
#   PFlog without cell-cycle-HVG filtering
#   PFlog with cell-cycle-HVG filtering
#
# Outputs:
#   Delegates to scripts/preprocess-sobj.R, which writes preprocessed Seurat
#   objects to CURRENT_OBJECT_DIR and preprocessing figures to FIGURE_DIR.

system("Rscript scripts/preprocess-sobj.R --normalization log1p")
system(
  "Rscript scripts/preprocess-sobj.R --normalization log1p --filter-cell-cycle"
)
system("Rscript scripts/preprocess-sobj.R --normalization pflog")
system(
  "Rscript scripts/preprocess-sobj.R --normalization pflog --filter-cell-cycle"
)
