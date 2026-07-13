#!/usr/bin/env Rscript

# Run every current preprocessing branch.
#
#   Rscript scripts/preprocess-all.R \
#     [--input <seurat-object.rds> | --input-source <legacy|counts-qc>]
#
# Arguments:
#   --input
#     Explicit Seurat object for every branch. Cannot be combined with
#     --input-source.
#   --input-source
#     Named source object for every branch. Defaults to legacy.
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

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) {
  i <- match(name, args)
  if (is.na(i)) {
    return(NULL)
  }
  if (i == length(args) || startsWith(args[[i + 1]], "--")) {
    return(TRUE)
  }
  args[[i + 1]]
}
arg_value <- function(name) {
  value <- arg(name)
  if (identical(value, TRUE)) {
    stop("Missing value for ", name, call. = FALSE)
  }
  value
}

input <- arg_value("--input")
input_source <- arg_value("--input-source")
if (!is.null(input) && !is.null(input_source)) {
  stop("Use either --input or --input-source, not both.", call. = FALSE)
}
input_args <- if (!is.null(input)) {
  c("--input", input)
} else if (!is.null(input_source)) {
  c("--input-source", input_source)
} else {
  character()
}
rscript <- file.path(R.home("bin"), "Rscript")

commands <- list(
  c(preprocess_script, input_args, "--normalization", "log1p"),
  c(
    preprocess_script,
    input_args,
    "--normalization",
    "log1p",
    "--filter-cell-cycle"
  ),
  c(preprocess_script, input_args, "--normalization", "pflog"),
  c(
    preprocess_script,
    input_args,
    "--normalization",
    "pflog",
    "--filter-cell-cycle"
  )
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
