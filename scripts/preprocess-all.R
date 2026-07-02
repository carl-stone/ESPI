system("Rscript scripts/preprocess-sobj.R --normalization log1p")
system(
  "Rscript scripts/preprocess-sobj.R --normalization log1p --filter-cell-cycle"
)
system("Rscript scripts/preprocess-sobj.R --normalization pflog")
system(
  "Rscript scripts/preprocess-sobj.R --normalization pflog --filter-cell-cycle"
)
