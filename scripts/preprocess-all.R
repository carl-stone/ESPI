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

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/preprocess-all.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

preprocess_script <- here::here("scripts", "preprocess-sobj.R")
rscript <- file.path(R.home("bin"), "Rscript")

commands <- list(
  c(preprocess_script, "--normalization", "log1p"),
  c(preprocess_script, "--normalization", "log1p", "--filter-cell-cycle"),
  c(preprocess_script, "--normalization", "pflog"),
  c(preprocess_script, "--normalization", "pflog", "--filter-cell-cycle")
)

# ---- work ----

for (command in commands) {
  status <- system2(rscript, command)
  if (!identical(as.integer(status), 0L)) {
    stop(
      "Preprocess command failed: ",
      paste(shQuote(c(rscript, command)), collapse = " "),
      call. = FALSE
    )
  }
}
