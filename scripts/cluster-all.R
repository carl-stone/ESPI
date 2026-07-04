#!/usr/bin/env Rscript

# Cluster every preprocessed branch in an input directory.
#
# Usage:
#   Rscript scripts/cluster-all.R \
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
#     Additional PC counts to cluster. Forwarded to cluster-sobj.R. Defaults
#     there to 30,50.
#   --resolutions
#     Leiden resolutions to cluster. Forwarded to cluster-sobj.R. Defaults
#     there to 0.3,0.5,0.8.
#   --dry-run
#     Print the Rscript commands that would run, without executing them.
#
# Outputs:
#   Delegates to scripts/cluster-sobj.R for every preprocess_*.rds input.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg) == 1) {
  dirname(dirname(normalizePath(
    sub("^--file=", "", file_arg),
    mustWork = TRUE
  )))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)

trailing_args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) {
  i <- match(name, trailing_args)
  if (is.na(i)) {
    return(NULL)
  }
  if (i == length(trailing_args) || startsWith(trailing_args[[i + 1]], "--")) {
    return(TRUE)
  }
  trailing_args[[i + 1]]
}
arg_flag <- function(name) {
  identical(arg(name), TRUE)
}

elbow_n <- arg("--elbow-n")
if (is.null(elbow_n)) {
  elbow_n <- "20"
}

suppressPackageStartupMessages({
  devtools::load_all(repo_root, export_all = FALSE, quiet = TRUE)
})

input_dir <- arg("--input-dir")
if (is.null(input_dir)) {
  input_dir <- CURRENT_OBJECT_DIR
}
input_dir <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)
if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir, call. = FALSE)
}

inputs <- list.files(
  input_dir,
  pattern = "^preprocess_.*\\.rds$",
  full.names = TRUE
)
if (length(inputs) == 0) {
  stop("No preprocessed objects found in: ", input_dir, call. = FALSE)
}

cluster_script <- file.path(repo_root, "scripts", "cluster-sobj.R")
rscript <- file.path(R.home("bin"), "Rscript")
dry_run <- arg_flag("--dry-run")

extra_dims <- arg("--extra-dims")
resolutions <- arg("--resolutions")

for (input in inputs) {
  command_args <- c(
    cluster_script,
    "--input",
    input,
    "--elbow-n",
    elbow_n
  )
  if (!is.null(extra_dims) && !identical(extra_dims, TRUE)) {
    command_args <- c(command_args, "--extra-dims", extra_dims)
  }
  if (!is.null(resolutions) && !identical(resolutions, TRUE)) {
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
