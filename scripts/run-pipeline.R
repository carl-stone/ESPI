#!/usr/bin/env Rscript

# Print the deterministic analysis plan or execute it through the figure stages.
#
# Usage:
#   Rscript scripts/run-pipeline.R --dry-run \
#     [--input-source counts-qc|legacy | --input <seurat.rds>] [--overwrite]
#
# Arguments:
#   --dry-run       Print the pipeline contract without checking or writing inputs.
#   --input-source  Use the named counts-qc or legacy source object.
#   --input         Use an explicit Seurat RDS input path.
#   --overwrite     Permit replacement of protected marker and DE outputs.

suppressPackageStartupMessages({
  library(here)
})
suppressMessages(here::i_am("scripts/run-pipeline.R"))
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

arguments <- commandArgs(trailingOnly = TRUE)
input_source <- "counts-qc"
input_path <- NULL
input_source_supplied <- FALSE
dry_run <- FALSE
overwrite <- FALSE

argument_index <- 1L
while (argument_index <= length(arguments)) {
  argument <- arguments[[argument_index]]

  if (identical(argument, "--dry-run")) {
    dry_run <- TRUE
    argument_index <- argument_index + 1L
    next
  }
  if (identical(argument, "--overwrite")) {
    overwrite <- TRUE
    argument_index <- argument_index + 1L
    next
  }
  if (argument %in% c("--input-source", "--input")) {
    if (
      argument_index == length(arguments) ||
        startsWith(arguments[[argument_index + 1L]], "--")
    ) {
      stop("Missing value for ", argument, ".", call. = FALSE)
    }

    value <- arguments[[argument_index + 1L]]
    if (identical(argument, "--input-source")) {
      input_source <- value
      input_source_supplied <- TRUE
    } else {
      input_path <- value
    }
    argument_index <- argument_index + 2L
    next
  }

  stop("Unknown argument: ", argument, ".", call. = FALSE)
}

if (!is.null(input_path) && input_source_supplied) {
  stop("Use either --input or --input-source, not both.", call. = FALSE)
}
if (!input_source %in% c("counts-qc", "legacy")) {
  stop(
    "--input-source must be one of counts-qc or legacy.",
    call. = FALSE
  )
}
if (!is.null(input_path)) {
  input_source <- "explicit"
} else {
  input_path <- if (identical(input_source, "counts-qc")) {
    file.path(INPUT_OBJECT_DIR, "sobj_qc_filtered.rds")
  } else {
    file.path(INPUT_OBJECT_DIR, "pipseq_processed_matrix_with_egfp.rds")
  }
}

run_spec <- list(
  input_source = input_source,
  input_path = input_path,
  overwrite = overwrite,
  normalization = "pflog",
  filter_cell_cycle_hvgs = FALSE,
  cluster_elbow_n = 20L,
  sensitivity_dims = c(30L, 50L),
  candidate_resolutions = c(0.3, 0.5, 0.8),
  source_chosen_dims = 30L,
  source_resolution = 0.3,
  mg_chosen_dims = 20L,
  mg_resolution = 0.5,
  mg_pca_dims = 50L,
  expression_layer = "pflog",
  marker_layer = "data",
  module_score_layer = "data",
  counts_layer = "counts"
)

# ---- pipeline plan ----

resolution_tag <- function(value) {
  format(value, trim = TRUE, scientific = FALSE)
}

cell_cycle_tag <- function(filtered) {
  if (isTRUE(filtered)) "filter_cc" else "no_filter_cc"
}

preprocess_path <- function(normalization, filtered) {
  file.path(
    CURRENT_OBJECT_DIR,
    sprintf(
      "preprocess_%s_%s.rds",
      normalization,
      if (isTRUE(filtered)) "filter-cc" else "no-filter-cc"
    )
  )
}

branch_tag <- function(normalization, filtered, dataset_tag = NULL) {
  paste(
    c(normalization, dataset_tag, cell_cycle_tag(filtered)),
    collapse = "_"
  )
}

cluster_path <- function(branch, elbow_n) {
  file.path(
    CURRENT_OBJECT_DIR,
    sprintf("cluster_%s_elbow%d.rds", branch, elbow_n)
  )
}

cluster_column <- function(branch, dims, resolution) {
  sprintf(
    "cluster_%s_dims%d_res%s",
    branch,
    dims,
    resolution_tag(resolution)
  )
}

cluster_expected_outputs <- function(branch) {
  candidate_dims <- c(run_spec$cluster_elbow_n, run_spec$sensitivity_dims)
  notebook_umap_paths <- as.vector(outer(
    candidate_dims,
    run_spec$candidate_resolutions,
    FUN = function(dims, resolution) {
      umap_name <- sprintf("umap_%s_dims%d", branch, dims)
      cluster_name <- cluster_column(branch, dims, resolution)
      file.path(
        here::here("notebook", "figures"),
        sprintf(
          "%s_by_%s.png",
          umap_name,
          gsub("[^A-Za-z0-9_-]", "_", cluster_name)
        )
      )
    }
  ))
  c(
    cluster_path(branch, run_spec$cluster_elbow_n),
    notebook_umap_paths
  )
}

heatmap_paths <- function(branch, type, output_dir, dims, resolution) {
  output_tag <- switch(
    type,
    marker = sprintf(
      "cell_type_marker_heatmap_%s_%s_cells_dims%d_res%s",
      run_spec$expression_layer,
      branch,
      dims,
      resolution_tag(resolution)
    ),
    module = sprintf(
      "cell_type_module_p27_heatmap_%s_%s_dims%d_res%s",
      run_spec$expression_layer,
      branch,
      dims,
      resolution_tag(resolution)
    )
  )
  png_path <- file.path(output_dir, paste0(output_tag, ".png"))
  pdf_path <- file.path(output_dir, paste0(output_tag, ".pdf"))
  notebook_png_path <- here::here(
    "notebook",
    "figures",
    paste0(output_tag, ".png")
  )
  if (identical(type, "marker")) {
    return(c(png_path, pdf_path, notebook_png_path))
  }
  table_dir <- file.path(TABLE_DIR, "annotation")
  c(
    png_path,
    pdf_path,
    file.path(table_dir, paste0(output_tag, "_module_scores.tsv")),
    file.path(table_dir, paste0(output_tag, "_p27_enrichment.tsv")),
    notebook_png_path
  )
}

mg_figure_paths <- function(branch, dims, resolution) {
  resolution_tag_value <- resolution_tag(resolution)
  output_dir <- file.path(FIGURE_DIR, "mg_selected")
  table_dir <- file.path(TABLE_DIR, "mg_selected")
  notebook_dir <- here::here("notebook", "figures")
  output_tags <- c(
    cluster = sprintf(
      "mg_selected_cluster_umap_%s_dims%d_res%s",
      branch,
      dims,
      resolution_tag_value
    ),
    condition = sprintf(
      "mg_selected_condition_umap_%s_dims%d_res%s",
      branch,
      dims,
      resolution_tag_value
    ),
    feature = sprintf(
      "mg_selected_feature_umap_%s_%s_dims%d_res%s",
      run_spec$expression_layer,
      branch,
      dims,
      resolution_tag_value
    ),
    coexpression = sprintf(
      "mg_selected_ascl1_hes6_coexpression_%s_dims%d_res%s",
      branch,
      dims,
      resolution_tag_value
    ),
    abundance = sprintf(
      "mg_selected_cluster_abundance_enrichment_%s_dims%d_res%s",
      branch,
      dims,
      resolution_tag_value
    ),
    proportion = sprintf(
      "mg_selected_cluster_proportion_by_mouse_%s_dims%d_res%s",
      branch,
      dims,
      resolution_tag_value
    )
  )
  c(
    unlist(
      lapply(
        output_tags,
        function(output_tag) {
          file.path(output_dir, paste0(output_tag, c(".png", ".pdf")))
        }
      ),
      use.names = FALSE
    ),
    file.path(
      table_dir,
      paste0(
        c(
          output_tags[["abundance"]],
          sprintf(
            "mg_selected_cluster_proportion_randomization_%s_dims%d_res%s",
            branch,
            dims,
            resolution_tag_value
          ),
          sprintf(
            "mg_selected_sample_cluster_proportions_%s_dims%d_res%s",
            branch,
            dims,
            resolution_tag_value
          )
        ),
        ".tsv"
      )
    ),
    file.path(notebook_dir, paste0(output_tags, ".png"))
  )
}

shell_quote <- function(argument) {
  paste0(
    "'",
    gsub("'", "'\"'\"'", as.character(argument), fixed = TRUE),
    "'"
  )
}

render_command <- function(command) {
  paste(vapply(command, shell_quote, character(1)), collapse = " ")
}

new_stage <- function(
  name,
  command,
  expects,
  protected_outputs = character(),
  validator = NULL
) {
  list(
    name = name,
    command = command,
    expects = expects,
    protected_outputs = protected_outputs,
    validator = validator
  )
}

read_stage_rds <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(error) {
      stop(
        "Cannot read expected RDS output ",
        path,
        ": ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
}

validate_seurat_object <- function(object, path) {
  if (!inherits(object, "Seurat")) {
    stop("Expected a Seurat object at ", path, ".", call. = FALSE)
  }
}

validate_counts_output <- function(stage) {
  raw_object <- read_stage_rds(stage$expects[["raw_counts"]])
  validate_seurat_object(raw_object, stage$expects[["raw_counts"]])

  missing_metadata <- setdiff(c("Mouse", "Condition"), colnames(raw_object[[]]))
  if (length(missing_metadata) > 0L) {
    stop(
      "Count output is missing required metadata column(s): ",
      paste(missing_metadata, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
}
validate_preprocess_outputs <- function(stage) {
  expected_branches <- names(stage$expects)
  if (
    is.null(expected_branches) ||
      any(!nzchar(expected_branches))
  ) {
    stop(
      "Preprocess output expectations must be named by branch.",
      call. = FALSE
    )
  }

  for (branch in expected_branches) {
    path <- stage$expects[[branch]]
    sobj <- read_stage_rds(path)
    validate_seurat_object(sobj, path)

    missing_metadata <- setdiff(c("Mouse", "Condition"), colnames(sobj[[]]))
    if (length(missing_metadata) > 0L) {
      stop(
        "Preprocess output ",
        path,
        " is missing required metadata column(s): ",
        paste(missing_metadata, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    pca_embeddings <- tryCatch(
      SeuratObject::Embeddings(sobj, reduction = "pca"),
      error = function(error) NULL
    )
    if (is.null(pca_embeddings) || ncol(pca_embeddings) != 50L) {
      stop(
        "Preprocess output ",
        path,
        " must contain a 50-PC PCA reduction.",
        call. = FALSE
      )
    }

    branch_parts <- strsplit(branch, "_", fixed = TRUE)[[1L]]
    expected_normalization <- branch_parts[[1L]]
    expected_filtered <- identical(branch_parts[[2L]], "filter")
    preprocessing <- sobj@misc$preprocessing
    if (
      !identical(preprocessing$normalization, expected_normalization) ||
        !identical(
          isTRUE(preprocessing$filtered_cell_cycle),
          expected_filtered
        )
    ) {
      stop(
        "Preprocess output ",
        path,
        " has branch metadata inconsistent with expected branch ",
        branch,
        ".",
        call. = FALSE
      )
    }
  }
}

validate_cluster_output <- function(
  stage,
  expected_branch,
  chosen_column = NULL
) {
  path <- stage$expects[[1L]]
  sobj <- read_stage_rds(path)
  validate_seurat_object(sobj, path)

  preprocessing <- sobj@misc$preprocessing
  branch_parts <- strsplit(expected_branch, "_", fixed = TRUE)[[1L]]
  expected_normalization <- branch_parts[[1L]]
  expected_filtered <- identical(branch_parts[[2L]], "filter")
  if (
    !identical(preprocessing$normalization, expected_normalization) ||
      !identical(
        isTRUE(preprocessing$filtered_cell_cycle),
        expected_filtered
      )
  ) {
    stop(
      "Clustered output ",
      path,
      " has preprocessing metadata inconsistent with branch ",
      expected_branch,
      ".",
      call. = FALSE
    )
  }

  clustering <- sobj@misc$clustering
  if (!identical(clustering$branch_tag, expected_branch)) {
    stop(
      "Clustered output ",
      path,
      " has branch tag ",
      paste(clustering$branch_tag, collapse = ", "),
      "; expected ",
      expected_branch,
      ".",
      call. = FALSE
    )
  }

  candidate_dims <- c(run_spec$cluster_elbow_n, run_spec$sensitivity_dims)
  candidate_columns <- as.vector(outer(
    candidate_dims,
    run_spec$candidate_resolutions,
    FUN = function(dims, resolution) {
      sprintf(
        "cluster_%s_dims%d_res%s",
        expected_branch,
        dims,
        resolution_tag(resolution)
      )
    }
  ))
  missing_candidates <- setdiff(candidate_columns, colnames(sobj[[]]))
  if (length(missing_candidates) > 0L) {
    stop(
      "Clustered output ",
      path,
      " is missing candidate cluster column(s): ",
      paste(missing_candidates, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  umap_reductions <- sprintf(
    "umap_%s_dims%d",
    expected_branch,
    candidate_dims
  )
  missing_umaps <- setdiff(umap_reductions, names(sobj@reductions))
  if (length(missing_umaps) > 0L) {
    stop(
      "Clustered output ",
      path,
      " is missing candidate UMAP reduction(s): ",
      paste(missing_umaps, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(chosen_column) && !chosen_column %in% colnames(sobj[[]])) {
    stop(
      "Chosen source output ",
      path,
      " is missing cluster column ",
      chosen_column,
      ".",
      call. = FALSE
    )
  }
}

make_cluster_validator <- function(expected_branch, chosen_column = NULL) {
  force(expected_branch)
  force(chosen_column)
  function(stage) {
    validate_cluster_output(stage, expected_branch, chosen_column)
  }
}


validate_mg_preprocess_outputs <- function(stage) {
  expected_branches <- names(mg_preprocess_paths)
  if (
    is.null(expected_branches) ||
      !identical(
        sort(expected_branches),
        sort(c("no_filter_cc", "filter_cc"))
      )
  ) {
    stop(
      "MG preprocess output expectations must contain both cell-cycle branches.",
      call. = FALSE
    )
  }

  expected_source_provenance <- list(
    source_cluster_column = run_spec$source_cluster_column,
    source_input = chosen_source_path,
    source_cluster_selection_table = mg_selection_table_path,
    source_cluster_selection_figure = mg_selection_diagnostic_paths[[1L]]
  )

  for (filter_key in expected_branches) {
    path <- mg_preprocess_paths[[filter_key]]
    sobj <- read_stage_rds(path)
    validate_seurat_object(sobj, path)

    missing_metadata <- setdiff(c("Mouse", "Condition"), colnames(sobj[[]]))
    if (length(missing_metadata) > 0L) {
      stop(
        "MG preprocess output ",
        path,
        " is missing required metadata column(s): ",
        paste(missing_metadata, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    pca_embeddings <- tryCatch(
      SeuratObject::Embeddings(sobj, reduction = "pca"),
      error = function(error) NULL
    )
    if (
      is.null(pca_embeddings) ||
        ncol(pca_embeddings) != run_spec$mg_pca_dims
    ) {
      stop(
        "MG preprocess output ",
        path,
        " must contain a ",
        run_spec$mg_pca_dims,
        "-PC PCA reduction.",
        call. = FALSE
      )
    }

    expected_filtered <- identical(filter_key, "filter_cc")
    preprocessing <- sobj@misc$preprocessing
    if (
      !identical(
        preprocessing$normalization,
        run_spec$normalization
      ) ||
        !identical(
          isTRUE(preprocessing$filtered_cell_cycle),
          expected_filtered
        ) ||
        !identical(
          preprocessing$dataset_tag,
          "mg_selected"
        )
    ) {
      stop(
        "MG preprocess output ",
        path,
        " has normalization, cell-cycle, or dataset metadata ",
        "inconsistent with branch ",
        filter_key,
        ".",
        call. = FALSE
      )
    }

    for (field in names(expected_source_provenance)) {
      if (
        !identical(
          preprocessing[[field]],
          expected_source_provenance[[field]]
        )
      ) {
        stop(
          "MG preprocess output ",
          path,
          " has source provenance field ",
          field,
          " inconsistent with the run specification.",
          call. = FALSE
        )
      }
    }
  }
}

validate_mg_cluster_output <- function(stage, expected_branch) {
  path <- stage$expects[[1L]]
  sobj <- read_stage_rds(path)
  validate_seurat_object(sobj, path)

  missing_metadata <- setdiff(c("Mouse", "Condition"), colnames(sobj[[]]))
  if (length(missing_metadata) > 0L) {
    stop(
      "MG clustered output ",
      path,
      " is missing required metadata column(s): ",
      paste(missing_metadata, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  preprocessing <- sobj@misc$preprocessing
  expected_filtered <- identical(
    expected_branch,
    branch_tag(
      run_spec$normalization,
      TRUE,
      dataset_tag = "mg_selected"
    )
  )
  if (
    !identical(preprocessing$normalization, run_spec$normalization) ||
      !identical(
        isTRUE(preprocessing$filtered_cell_cycle),
        expected_filtered
      ) ||
      !identical(preprocessing$dataset_tag, "mg_selected")
  ) {
    stop(
      "MG clustered output ",
      path,
      " has preprocessing metadata inconsistent with branch ",
      expected_branch,
      ".",
      call. = FALSE
    )
  }

  expected_source_provenance <- list(
    source_cluster_column = run_spec$source_cluster_column,
    source_input = chosen_source_path,
    source_cluster_selection_table = mg_selection_table_path,
    source_cluster_selection_figure = mg_selection_diagnostic_paths[[1L]]
  )
  for (field in names(expected_source_provenance)) {
    if (
      !identical(
        preprocessing[[field]],
        expected_source_provenance[[field]]
      )
    ) {
      stop(
        "MG clustered output ",
        path,
        " has source provenance field ",
        field,
        " inconsistent with the run specification.",
        call. = FALSE
      )
    }
  }

  clustering <- sobj@misc$clustering
  if (!identical(clustering$branch_tag, expected_branch)) {
    stop(
      "MG clustered output ",
      path,
      " has branch tag ",
      paste(clustering$branch_tag, collapse = ", "),
      "; expected ",
      expected_branch,
      ".",
      call. = FALSE
    )
  }

  candidate_dims <- c(run_spec$cluster_elbow_n, run_spec$sensitivity_dims)
  candidate_columns <- as.vector(outer(
    candidate_dims,
    run_spec$candidate_resolutions,
    FUN = function(dims, resolution) {
      sprintf(
        "cluster_%s_dims%d_res%s",
        expected_branch,
        dims,
        resolution_tag(resolution)
      )
    }
  ))
  missing_candidates <- setdiff(candidate_columns, colnames(sobj[[]]))
  if (length(missing_candidates) > 0L) {
    stop(
      "MG clustered output ",
      path,
      " is missing candidate cluster column(s): ",
      paste(missing_candidates, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  declared_candidates <- clustering$candidate_names
  if (
    !is.character(declared_candidates) ||
      length(setdiff(candidate_columns, declared_candidates)) > 0L
  ) {
    stop(
      "MG clustered output ",
      path,
      " does not declare all candidate cluster columns.",
      call. = FALSE
    )
  }

  umap_reductions <- sprintf(
    "umap_%s_dims%d",
    expected_branch,
    candidate_dims
  )
  missing_umaps <- setdiff(umap_reductions, names(sobj@reductions))
  if (length(missing_umaps) > 0L) {
    stop(
      "MG clustered output ",
      path,
      " is missing candidate UMAP reduction(s): ",
      paste(missing_umaps, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  expected_chosen_column <- sprintf(
    "cluster_%s_dims%d_res%s",
    expected_branch,
    run_spec$mg_chosen_dims,
    resolution_tag(run_spec$mg_resolution)
  )
  if (!expected_chosen_column %in% colnames(sobj[[]])) {
    stop(
      "MG clustered output ",
      path,
      " is missing chosen cluster column ",
      expected_chosen_column,
      ".",
      call. = FALSE
    )
  }
}

make_mg_cluster_validator <- function(expected_branch) {
  force(expected_branch)
  function(stage) validate_mg_cluster_output(stage, expected_branch)
}

validate_qc_outputs <- function(stage) {
  annotated_path <- stage$expects[["annotated_raw"]]
  filtered_path <- stage$expects[["filtered"]]
  annotated_object <- read_stage_rds(annotated_path)
  filtered_object <- read_stage_rds(filtered_path)
  validate_seurat_object(annotated_object, annotated_path)
  validate_seurat_object(filtered_object, filtered_path)

  annotated_metadata <- annotated_object[[]]
  filtered_metadata <- filtered_object[[]]
  missing_annotation_columns <- setdiff(
    c("is_cell", "pass_qc"),
    colnames(annotated_metadata)
  )
  if (length(missing_annotation_columns) > 0L) {
    stop(
      "Annotated QC output is missing required metadata column(s): ",
      paste(missing_annotation_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.logical(annotated_metadata$pass_qc)) {
    stop("Annotated QC output pass_qc metadata must be logical.", call. = FALSE)
  }

  missing_filtered_columns <- setdiff(
    c("Mouse", "Condition"),
    colnames(filtered_metadata)
  )
  if (length(missing_filtered_columns) > 0L) {
    stop(
      "Filtered QC output is missing required metadata column(s): ",
      paste(missing_filtered_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  expected_cells <- colnames(annotated_object)[
    annotated_metadata$pass_qc %in% TRUE
  ]
  if (!identical(colnames(filtered_object), expected_cells)) {
    stop(
      "Filtered QC cells do not exactly match annotated cells with pass_qc TRUE.",
      call. = FALSE
    )
  }
}

validate_stage_outputs <- function(stage) {
  missing_paths <- stage$expects[!file.exists(stage$expects)]
  if (length(missing_paths) > 0L) {
    stop(
      "Stage ",
      stage$name,
      " did not create expected output(s): ",
      paste(missing_paths, collapse = "; "),
      call. = FALSE
    )
  }
  if (is.function(stage$validator)) {
    stage$validator(stage)
  }
}

run_stage <- function(stage) {
  message("Running stage ", stage$name, ": ", render_command(stage$command))
  status <- tryCatch(
    suppressWarnings(
      system2(
        command = stage$command[[1]],
        args = stage$command[-1],
        stdout = "",
        stderr = ""
      )
    ),
    error = function(error) {
      stop(
        "Stage ",
        stage$name,
        " could not start: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (!identical(status, 0L)) {
    stop(
      "Stage ",
      stage$name,
      " failed with exit status ",
      status,
      ".",
      call. = FALSE
    )
  }
  tryCatch(
    validate_stage_outputs(stage),
    error = function(error) {
      stop(
        "Stage ",
        stage$name,
        " output validation failed: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  message("Completed stage ", stage$name, ".")
}

source_branches <- data.frame(
  normalization = c("log1p", "log1p", "pflog", "pflog"),
  filtered = c(FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
source_branches$branch <- vapply(
  seq_len(nrow(source_branches)),
  function(index) {
    branch_tag(
      source_branches$normalization[[index]],
      source_branches$filtered[[index]]
    )
  },
  character(1)
)
source_branches$preprocess_path <- vapply(
  seq_len(nrow(source_branches)),
  function(index) {
    preprocess_path(
      source_branches$normalization[[index]],
      source_branches$filtered[[index]]
    )
  },
  character(1)
)

run_spec$source_branch_tag <- branch_tag(
  run_spec$normalization,
  run_spec$filter_cell_cycle_hvgs
)
run_spec$mg_branch_tag <- branch_tag(
  run_spec$normalization,
  run_spec$filter_cell_cycle_hvgs,
  dataset_tag = "mg_selected"
)
run_spec$source_cluster_column <- cluster_column(
  run_spec$source_branch_tag,
  run_spec$source_chosen_dims,
  run_spec$source_resolution
)
run_spec$mg_cluster_column <- cluster_column(
  run_spec$mg_branch_tag,
  run_spec$mg_chosen_dims,
  run_spec$mg_resolution
)

chosen_source_path <- cluster_path(
  run_spec$source_branch_tag,
  run_spec$cluster_elbow_n
)
mg_preprocess_paths <- c(
  no_filter_cc = file.path(
    CURRENT_OBJECT_DIR,
    sprintf(
      "preprocess_%s_mg_selected_no-filter-cc.rds",
      run_spec$normalization
    )
  ),
  filter_cc = file.path(
    CURRENT_OBJECT_DIR,
    sprintf(
      "preprocess_%s_mg_selected_filter-cc.rds",
      run_spec$normalization
    )
  )
)
mg_selection_table_path <- file.path(
  TABLE_DIR,
  "mg_selected",
  "mg_selected_cluster_selection.tsv"
)
mg_selection_diagnostic_paths <- file.path(
  FIGURE_DIR,
  "mg_selected",
  paste0("mg_selected_cluster_selection_diagnostics", c(".png", ".pdf"))
)
mg_elbow_expected_outputs <- unlist(
  lapply(
    c("no_filter_cc", "filter_cc"),
    function(filter_key) {
      file.path(
        FIGURE_DIR,
        "mg_selected",
        paste0(
          "elbow_pflog_mg_selected_",
          filter_key,
          c(".png", ".pdf")
        )
      )
    }
  ),
  use.names = FALSE
)
mg_selection_expected_outputs <- c(
  unname(mg_preprocess_paths),
  mg_selection_table_path,
  mg_selection_diagnostic_paths,
  mg_elbow_expected_outputs
)
mg_branch_paths <- setNames(
  vapply(
    c(FALSE, TRUE),
    function(filtered) {
      cluster_path(
        branch_tag(
          run_spec$normalization,
          filtered,
          dataset_tag = "mg_selected"
        ),
        run_spec$cluster_elbow_n
      )
    },
    character(1)
  ),
  c("no_filter_cc", "filter_cc")
)
mg_summary_expected_outputs <- c(
  file.path(
    TABLE_DIR,
    "mg_selected",
    "mg_selected_cluster_grid_summary.tsv"
  ),
  unlist(
    lapply(
      c("no_filter_cc", "filter_cc"),
      function(filter_key) {
        file.path(
          FIGURE_DIR,
          "mg_selected",
          paste0(
            "mg_selected_umap_resolution_sweep_pflog_mg_selected_",
            filter_key,
            "_dims",
            run_spec$mg_pca_dims,
            c(".png", ".pdf")
          )
        )
      }
    ),
    use.names = FALSE
  )
)

marker_table_dir <- file.path(TABLE_DIR, "mg_selected")
marker_figure_dir <- file.path(FIGURE_DIR, "mg_selected")
marker_tag <- sprintf(
  "data_%s_dims%d_res%s",
  run_spec$mg_branch_tag,
  run_spec$mg_chosen_dims,
  resolution_tag(run_spec$mg_resolution)
)
marker_protected_outputs <- c(
  file.path(marker_table_dir, paste0("find_all_markers_", marker_tag, ".csv")),
  file.path(
    marker_table_dir,
    sprintf(
      "find_all_markers_top5_%s.csv",
      marker_tag
    )
  ),
  file.path(
    marker_table_dir,
    sprintf(
      "find_all_markers_summary_%s.csv",
      marker_tag
    )
  ),
  file.path(
    marker_table_dir,
    sprintf(
      "find_all_markers_identity_map_%s_dims%d_res%s.csv",
      run_spec$mg_branch_tag,
      run_spec$mg_chosen_dims,
      resolution_tag(run_spec$mg_resolution)
    )
  ),
  file.path(
    marker_figure_dir,
    sprintf("mg_selected_cluster_marker_dotplot_%s_top5.png", marker_tag)
  ),
  file.path(
    marker_figure_dir,
    sprintf("mg_selected_cluster_marker_dotplot_%s_top5.pdf", marker_tag)
  ),
  file.path(
    here::here("notebook", "figures"),
    sprintf(
      "mg_selected_cluster_marker_dotplot_%s_top5.png",
      marker_tag
    )
  )
)

de_output_dir <- file.path(DEG_DIR, "mg_selected")
enrichment_output_dir <- file.path(ENRICHMENT_DIR, "mg_selected")
de_figure_dir <- file.path(FIGURE_DIR, "mg_selected")
de_protected_outputs <- c(
  file.path(
    de_output_dir,
    c(
      "pseudobulk_sample_summary.tsv",
      "design_summary.tsv",
      "deseq2_full_results.tsv",
      "deseq2_significant_degs.tsv",
      "deseq2_marker_overlap.tsv",
      "deseq2_paired_sensitivity_full_results.tsv",
      "deseq2_paired_sensitivity_significant_degs.tsv",
      "deseq2_paired_sensitivity_marker_overlap.tsv",
      "numbers.json"
    )
  ),
  file.path(
    enrichment_output_dir,
    c(
      "go_bp_ora_up.tsv",
      "go_bp_ora_down.tsv",
      "go_bp_gsea.tsv",
      "go_bp_gsea_symbol_entrez_mapping.tsv",
      "go_bp_ora_up_simplified.tsv",
      "go_bp_ora_down_simplified.tsv",
      "go_bp_gsea_simplified.tsv",
      "go_bp_ora_up_bayes_simplified.tsv",
      "go_bp_ora_down_bayes_simplified.tsv"
    )
  ),
  file.path(
    de_figure_dir,
    c(
      "mg_selected_de_volcano.png",
      "mg_selected_de_volcano.pdf",
      "mg_selected_go_ora_up_dotplot.png",
      "mg_selected_go_ora_up_dotplot.pdf",
      "mg_selected_go_ora_down_dotplot.png",
      "mg_selected_go_ora_down_dotplot.pdf",
      "mg_selected_go_gsea_dotplot.png",
      "mg_selected_go_gsea_dotplot.pdf",
      "mg_selected_go_ora_up_bayes_dotplot.png",
      "mg_selected_go_ora_up_bayes_dotplot.pdf",
      "mg_selected_go_ora_down_bayes_dotplot.png",
      "mg_selected_go_ora_down_bayes_dotplot.pdf"
    )
  ),
  file.path(
    here::here("notebook", "figures"),
    c(
      "mg_selected_de_volcano.png",
      "mg_selected_go_ora_up_dotplot.png",
      "mg_selected_go_ora_down_dotplot.png",
      "mg_selected_go_gsea_dotplot.png",
      "mg_selected_go_ora_up_bayes_dotplot.png",
      "mg_selected_go_ora_down_bayes_dotplot.png"
    )
  )
)

source_input_command <- if (identical(run_spec$input_source, "explicit")) {
  c("--input", run_spec$input_path)
} else {
  c("--input-source", run_spec$input_source)
}
overwrite_argument <- if (isTRUE(run_spec$overwrite)) {
  "--overwrite"
} else {
  character()
}

source_summary_expected_outputs <- c(
  file.path(
    TABLE_DIR,
    "cluster",
    c(
      "cluster_grid_summary.tsv",
      "cluster_grid_stability_summary.tsv",
      "cluster_grid_pairwise_stability.tsv"
    )
  ),
  file.path(
    FIGURE_DIR,
    "cluster",
    c(
      "cluster_grid_clustree_12_panel.png",
      "cluster_grid_clustree_12_panel.pdf",
      "umap_resolution_sweep_pflog_filter_cc_dims50.png",
      "umap_resolution_sweep_pflog_filter_cc_dims50.pdf"
    )
  )
)


stage_plan <- list()
if (identical(run_spec$input_source, "counts-qc")) {
  raw_counts_dir <- file.path(DATA_ROOT_DIR, "data", "input", "Raw Matrices")
  stage_plan <- c(
    stage_plan,
    list(
      new_stage(
        "process-counts",
        c("Rscript", "scripts/01-process-counts.R"),
        c(
          raw_counts = file.path(
            DATA_ROOT_DIR,
            "data",
            "input",
            "sobj_raw.rds"
          )
        ),
        validator = validate_counts_output
      ),
      new_stage(
        "qc-filtering",
        c("Rscript", "scripts/02-qc-filtering.R"),
        c(
          annotated_raw = file.path(
            INPUT_OBJECT_DIR,
            "sobj_raw_with_qc.rds"
          ),
          filtered = file.path(INPUT_OBJECT_DIR, "sobj_qc_filtered.rds")
        ),
        validator = validate_qc_outputs
      )
    )
  )
}
stage_plan <- c(
  stage_plan,
  list(
    new_stage(
      "preprocess-source",
      c(
        "Rscript",
        "scripts/03-preprocess-all.R",
        source_input_command
      ),
      setNames(
        unname(source_branches$preprocess_path),
        source_branches$branch
      ),
      validator = validate_preprocess_outputs
    )
  )
)

for (index in seq_len(nrow(source_branches))) {
  source_branch <- source_branches[index, ]
  stage_plan <- c(
    stage_plan,
    list(
      new_stage(
        sprintf(
          "cluster-source-%s-%s",
          source_branch$normalization,
          if (source_branch$filtered) "filter-cc" else "no-filter-cc"
        ),
        c(
          "Rscript",
          "scripts/04-cluster.R",
          "--input",
          source_branch$preprocess_path,
          "--elbow-n",
          as.character(run_spec$cluster_elbow_n),
          "--extra-dims",
          paste(run_spec$sensitivity_dims, collapse = ","),
          "--resolutions",
          paste(resolution_tag(run_spec$candidate_resolutions), collapse = ",")
        ),
        cluster_expected_outputs(source_branch$branch),
        validator = make_cluster_validator(
          source_branch$branch,
          if (identical(source_branch$branch, run_spec$source_branch_tag)) {
            run_spec$source_cluster_column
          } else {
            NULL
          }
        )
      )
    )
  )
}

stage_plan <- c(
  stage_plan,
  list(
    new_stage(
      "summarize-source",
      c("Rscript", "scripts/05-summarize-clusters.R"),
      source_summary_expected_outputs
    ),
    new_stage(
      "select-mg",
      c(
        "Rscript",
        "scripts/07-select-mg-subset.R",
        "--input",
        chosen_source_path,
        "--cluster-column",
        run_spec$source_cluster_column,
        "--dims",
        as.character(run_spec$mg_pca_dims),
        "--score-layer",
        run_spec$module_score_layer,
        "--dataset-tag",
        "mg_selected"
      ),
      mg_selection_expected_outputs,
      validator = validate_mg_preprocess_outputs
    )
  )
)

for (filter_key in names(mg_preprocess_paths)) {
  filtered <- identical(filter_key, "filter_cc")
  mg_branch <- branch_tag(
    run_spec$normalization,
    filtered,
    dataset_tag = "mg_selected"
  )
  stage_plan <- c(
    stage_plan,
    list(
      new_stage(
        paste0("cluster-mg-", if (filtered) "filter-cc" else "no-filter-cc"),
        c(
          "Rscript",
          "scripts/04-cluster.R",
          "--input",
          mg_preprocess_paths[[filter_key]],
          "--elbow-n",
          as.character(run_spec$cluster_elbow_n),
          "--extra-dims",
          paste(run_spec$sensitivity_dims, collapse = ","),
          "--resolutions",
          paste(resolution_tag(run_spec$candidate_resolutions), collapse = ",")
        ),
        cluster_expected_outputs(mg_branch),
        validator = make_mg_cluster_validator(mg_branch)
      )
    )
  )
}
stage_plan <- c(
  stage_plan,
  list(
    new_stage(
      "summarize-mg",
      c(
        "Rscript",
        "scripts/08-summarize-mg-clusters.R",
        "--elbow-n",
        as.character(run_spec$cluster_elbow_n)
      ),
      mg_summary_expected_outputs
    )
  )
)

source_marker_paths <- heatmap_paths(
  run_spec$source_branch_tag,
  "marker",
  file.path(FIGURE_DIR, "annotation"),
  run_spec$source_chosen_dims,
  run_spec$source_resolution
)
source_module_paths <- heatmap_paths(
  run_spec$source_branch_tag,
  "module",
  file.path(FIGURE_DIR, "annotation"),
  run_spec$source_chosen_dims,
  run_spec$source_resolution
)
stage_plan <- c(
  stage_plan,
  list(
    new_stage(
      "marker-heatmap-source",
      c(
        "Rscript",
        "scripts/06-plot-marker-heatmap.R",
        "--input",
        chosen_source_path,
        "--dims",
        as.character(run_spec$source_chosen_dims),
        "--resolution",
        resolution_tag(run_spec$source_resolution),
        "--layer",
        run_spec$expression_layer,
        "--out-dir",
        file.path(FIGURE_DIR, "annotation")
      ),
      source_marker_paths
    )
  )
)
for (filter_key in names(mg_branch_paths)) {
  filtered <- identical(filter_key, "filter_cc")
  mg_branch <- branch_tag(
    run_spec$normalization,
    filtered,
    dataset_tag = "mg_selected"
  )
  marker_paths <- heatmap_paths(
    mg_branch,
    "marker",
    file.path(FIGURE_DIR, "annotation"),
    run_spec$mg_chosen_dims,
    run_spec$mg_resolution
  )
  stage_plan <- c(
    stage_plan,
    list(
      new_stage(
        paste0(
          "marker-heatmap-mg-",
          if (filtered) "filter-cc" else "no-filter-cc"
        ),
        c(
          "Rscript",
          "scripts/06-plot-marker-heatmap.R",
          "--input",
          mg_branch_paths[[filter_key]],
          "--dims",
          as.character(run_spec$mg_chosen_dims),
          "--resolution",
          resolution_tag(run_spec$mg_resolution),
          "--layer",
          run_spec$expression_layer,
          "--out-dir",
          file.path(FIGURE_DIR, "annotation")
        ),
        marker_paths
      )
    )
  )
}

stage_plan <- c(
  stage_plan,
  list(
    new_stage(
      "module-heatmap-source",
      c(
        "Rscript",
        "scripts/10-plot-cluster-marker-heatmaps.R",
        "--input",
        chosen_source_path,
        "--dims",
        as.character(run_spec$source_chosen_dims),
        "--resolution",
        resolution_tag(run_spec$source_resolution),
        "--layer",
        run_spec$expression_layer,
        "--slot",
        run_spec$module_score_layer,
        "--out-dir",
        file.path(FIGURE_DIR, "annotation")
      ),
      source_module_paths
    )
  )
)
for (filter_key in names(mg_branch_paths)) {
  filtered <- identical(filter_key, "filter_cc")
  mg_branch <- branch_tag(
    run_spec$normalization,
    filtered,
    dataset_tag = "mg_selected"
  )
  module_paths <- heatmap_paths(
    mg_branch,
    "module",
    file.path(FIGURE_DIR, "annotation"),
    run_spec$mg_chosen_dims,
    run_spec$mg_resolution
  )
  stage_plan <- c(
    stage_plan,
    list(
      new_stage(
        paste0(
          "module-heatmap-mg-",
          if (filtered) "filter-cc" else "no-filter-cc"
        ),
        c(
          "Rscript",
          "scripts/10-plot-cluster-marker-heatmaps.R",
          "--input",
          mg_branch_paths[[filter_key]],
          "--dims",
          as.character(run_spec$mg_chosen_dims),
          "--resolution",
          resolution_tag(run_spec$mg_resolution),
          "--layer",
          run_spec$expression_layer,
          "--slot",
          run_spec$module_score_layer,
          "--out-dir",
          file.path(FIGURE_DIR, "annotation")
        ),
        module_paths
      )
    )
  )
}

for (filter_key in names(mg_branch_paths)) {
  filtered <- identical(filter_key, "filter_cc")
  mg_branch <- branch_tag(
    run_spec$normalization,
    filtered,
    dataset_tag = "mg_selected"
  )
  stage_plan <- c(
    stage_plan,
    list(
      new_stage(
        paste0("mg-figures-", if (filtered) "filter-cc" else "no-filter-cc"),
        c(
          "Rscript",
          "scripts/09-plot-mg-figures.R",
          "--input",
          mg_branch_paths[[filter_key]],
          "--branch-tag",
          mg_branch,
          "--elbow-n",
          as.character(run_spec$cluster_elbow_n),
          "--dims",
          as.character(run_spec$mg_chosen_dims),
          "--resolution",
          resolution_tag(run_spec$mg_resolution),
          "--layer",
          run_spec$expression_layer
        ),
        mg_figure_paths(
          mg_branch,
          run_spec$mg_chosen_dims,
          run_spec$mg_resolution
        )
      )
    )
  )
}

stage_plan <- c(
  stage_plan,
  list(
    new_stage(
      "mg-markers",
      c(
        "Rscript",
        "scripts/11-find-mg-markers.R",
        "--input",
        mg_branch_paths[["no_filter_cc"]],
        "--branch-tag",
        run_spec$mg_branch_tag,
        "--elbow-n",
        as.character(run_spec$cluster_elbow_n),
        "--dims",
        as.character(run_spec$mg_chosen_dims),
        "--resolution",
        resolution_tag(run_spec$mg_resolution),
        "--layer",
        run_spec$marker_layer,
        "--counts-layer",
        run_spec$counts_layer,
        "--confirm-no-merge",
        overwrite_argument
      ),
      marker_protected_outputs,
      protected_outputs = marker_protected_outputs
    ),
    new_stage(
      "mg-de",
      c(
        "Rscript",
        "scripts/12-run-mg-de.R",
        "--input",
        mg_branch_paths[["no_filter_cc"]],
        "--cluster-column",
        run_spec$mg_cluster_column,
        "--counts-layer",
        run_spec$counts_layer,
        "--lfc-shrink-type",
        "apeglm",
        overwrite_argument
      ),
      de_protected_outputs,
      protected_outputs = de_protected_outputs
    ),
    new_stage(
      "render-notebook",
      c("quarto", "render", "notebook/sc_analysis.qmd"),
      here::here("notebook", "sc_analysis.html")
    ),
    new_stage(
      "tripwires",
      c("Rscript", "tools/run-tripwires.R"),
      here::here("tools", "run-tripwires.R")
    )
  )
)

# ---- validation ----

preflight_source_inputs <- function() {
  required_paths <- if (identical(run_spec$input_source, "counts-qc")) {
    raw_counts_dir <- file.path(DATA_ROOT_DIR, "data", "input", "Raw Matrices")
    c(raw_counts_dir, file.path(raw_counts_dir, "Sample_Metadata_MS1.txt"))
  } else {
    run_spec$input_path
  }
  exists <- if (identical(run_spec$input_source, "counts-qc")) {
    c(
      dir.exists(required_paths[[1]]),
      file.exists(required_paths[[2]])
    )
  } else {
    file.exists(required_paths)
  }
  missing_paths <- required_paths[!exists]
  if (length(missing_paths) > 0L) {
    stop(
      "Pipeline input(s) do not exist: ",
      paste(missing_paths, collapse = "; "),
      call. = FALSE
    )
  }
}

preflight_protected_outputs <- function(stages) {
  protected_paths <- unlist(
    lapply(stages, `[[`, "protected_outputs"),
    use.names = FALSE
  )
  existing_paths <- protected_paths[
    file.exists(protected_paths) | nzchar(Sys.readlink(protected_paths))
  ]
  if (length(existing_paths) > 0L && !isTRUE(run_spec$overwrite)) {
    stop(
      "Protected pipeline output(s) already exist; use --overwrite to replace: ",
      paste(existing_paths, collapse = "; "),
      call. = FALSE
    )
  }
}

# ---- dry run ----

if (dry_run) {
  contract_lines <- c(
    "mode: dry-run",
    paste0("input_source: ", run_spec$input_source),
    paste0("overwrite: ", tolower(as.character(run_spec$overwrite))),
    paste0("source_cluster_column: ", run_spec$source_cluster_column),
    paste0("mg_cluster_column: ", run_spec$mg_cluster_column),
    paste0("mg_pca_dims: ", run_spec$mg_pca_dims),
    paste0("first_stage: ", stage_plan[[1]]$name),
    "final_stage: tripwires",
    paste0("stage_count: ", length(stage_plan))
  )
  base::cat(contract_lines, sep = "\n")
  for (stage in stage_plan) {
    base::cat(
      paste0("stage: ", stage$name),
      paste0("command: ", render_command(stage$command)),
      paste0("expects: ", paste(stage$expects, collapse = ";")),
      sep = "\n"
    )
  }
  quit(status = 0L, save = "no")
}

source_execution_stage_names <- c(
  "preprocess-source",
  paste0(
    "cluster-source-",
    source_branches$normalization,
    "-",
    ifelse(source_branches$filtered, "filter-cc", "no-filter-cc")
  ),
  "summarize-source"
)
mg_execution_stage_names <- c(
  "select-mg",
  "cluster-mg-no-filter-cc",
  "cluster-mg-filter-cc",
  "summarize-mg"
)
figure_execution_stage_names <- c(
  "marker-heatmap-source",
  "marker-heatmap-mg-no-filter-cc",
  "marker-heatmap-mg-filter-cc",
  "module-heatmap-source",
  "module-heatmap-mg-no-filter-cc",
  "module-heatmap-mg-filter-cc",
  "mg-figures-no-filter-cc",
  "mg-figures-filter-cc"
)
analysis_execution_stage_names <- c(
  "mg-markers",
  "mg-de",
  "render-notebook",
  "tripwires"
)
execution_stage_names <- if (identical(run_spec$input_source, "counts-qc")) {
  c(
    "process-counts",
    "qc-filtering",
    source_execution_stage_names,
    mg_execution_stage_names,
    figure_execution_stage_names,
    analysis_execution_stage_names
  )
} else {
  c(
    source_execution_stage_names,
    mg_execution_stage_names,
    figure_execution_stage_names,
    analysis_execution_stage_names
  )
}
execution_stages <- stage_plan[
  match(execution_stage_names, vapply(stage_plan, `[[`, character(1), "name"))
]
if (any(vapply(execution_stages, is.null, logical(1)))) {
  stop(
    "Source execution stages are not available in the pipeline plan.",
    call. = FALSE
  )
}

preflight_source_inputs()
preflight_protected_outputs(execution_stages)
for (stage in execution_stages) {
  run_stage(stage)
}
message("Pipeline execution completed successfully through tripwires.")
quit(status = 0L, save = "no")
