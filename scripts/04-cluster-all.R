#!/usr/bin/env Rscript

# Cluster every preprocessed branch in an input directory.
#
# Usage:
#   Rscript scripts/04-cluster-all.R \
#     --elbow-n <positive integer> \
#     --input-dir <directory> \
#     --extra-dims <comma-separated integers> \
#     --resolutions <comma-separated numbers> \
#     --dry-run
#
# Arguments:
#   --elbow-n
#     Primary PC count selected from elbow diagnostics. Defaults to 20.
#   --input-dir
#     Directory containing preprocess_*.rds inputs. Defaults to
#     CURRENT_OBJECT_DIR.
#   --extra-dims
#     Additional PC counts to cluster. Forwarded to 04-cluster.R. Defaults
#     there to 30,50.
#   --resolutions
#     Leiden resolutions to cluster. Forwarded to 04-cluster.R. Defaults
#     there to 0.3,0.5,0.8.
#   --dry-run
#     Print the Rscript commands that would run, without executing them.
#
# Outputs:
#   Delegates to scripts/04-cluster.R for every preprocess_*.rds input.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/04-cluster-all.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

args <- commandArgs(trailingOnly = TRUE)
# ANALYSIS_OK[R025]: wrapper keeps a local commandArgs() parser for
# independently runnable RStudio steps; CLI tripwires exercise the contract.
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
# ANALYSIS_OK[R025]: local value parser is intentionally duplicated for
# standalone wrapper execution; CLI tripwires exercise this narrow contract.
arg_value <- function(name, default = NULL, required = FALSE) {
  value <- arg(name)
  if (identical(value, TRUE)) {
    stop("Missing value for ", name, call. = FALSE)
  }
  if (is.null(value)) {
    if (required) {
      stop("Missing required argument ", name, call. = FALSE)
    }
    return(default)
  }
  value
}
# ANALYSIS_OK[R025]: local boolean parser remains self-contained for
# RStudio-step-friendly execution; CLI tripwires exercise the contract.
arg_flag <- function(name) {
  identical(arg(name), TRUE)
}

elbow_n <- arg_value("--elbow-n", default = "20")
input_dir <- arg_value("--input-dir", default = CURRENT_OBJECT_DIR)
extra_dims <- arg_value("--extra-dims", default = NULL)
resolutions <- arg_value("--resolutions", default = NULL)
dry_run <- arg_flag("--dry-run")

input_dir <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)
cluster_script <- here::here("scripts", "04-cluster.R")
rscript <- file.path(R.home("bin"), "Rscript")

# ---- validation ----

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir, call. = FALSE)
}

# ---- work ----

inputs <- list.files(
  input_dir,
  pattern = "^preprocess_.*\\.rds$",
  full.names = TRUE
)
if (length(inputs) == 0) {
  stop("No preprocessed objects found in: ", input_dir, call. = FALSE)
}

# ---- output ----

for (input in inputs) {
  command_args <- c(
    cluster_script,
    "--input",
    input,
    "--elbow-n",
    elbow_n
  )
  if (!is.null(extra_dims)) {
    command_args <- c(command_args, "--extra-dims", extra_dims)
  }
  if (!is.null(resolutions)) {
    command_args <- c(command_args, "--resolutions", resolutions)
  }

  if (dry_run) {
    message(paste(shQuote(c(rscript, command_args)), collapse = " "))
  } else {
    status <- system2(rscript, command_args)
    if (!identical(as.integer(status), 0L)) {
      stop(
        "Cluster command failed for input: ",
        input,
        " (exit status ",
        status,
        ")",
        call. = FALSE
      )
    }
  }
}
