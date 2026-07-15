#!/usr/bin/env Rscript

# ESPI tripwire runner.
# All checks use base R and avoid running the full analysis or rendering Quarto.

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
make_result <- function(status, slug, message) {
  list(status = status, slug = slug, message = message)
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
pass <- function(slug, message) make_result("PASS", slug, message)
# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
fail <- function(slug, message) make_result("FAIL", slug, message)
# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
skip <- function(slug, message) make_result("SKIP", slug, message)

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
script_path <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  sub("^--file=", "", file_arg[[1]])
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
find_repo_root <- function() {
  markers <- c(
    "analysis_labels.yml",
    file.path("scripts", "04-cluster.R"),
    file.path("notebook", "sc_analysis.qmd")
  )
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (all(file.exists(file.path(cwd, markers)))) {
    return(cwd)
  }

  sp <- script_path()
  if (!is.na(sp)) {
    sp <- normalizePath(sp, winslash = "/", mustWork = TRUE)
    candidate <- dirname(dirname(sp))
    if (all(file.exists(file.path(candidate, markers)))) {
      return(candidate)
    }
  }

  NA_character_
}

read_text <- function(path) {
  readLines(path, warn = FALSE)
}

squash <- function(x) {
  paste(x, collapse = " ")
}

# ANALYSIS_OK[script-entrypoint]: helper is an internal entrypoint exercised by main(); R026 cross-file caller checks do not model this script's dispatch graph.
compact_problem_list <- function(problems, max_n = MAX_PROBLEM_SUMMARY_ITEMS) {
  if (length(problems) <= max_n) {
    return(paste(problems, collapse = " | "))
  }
  paste(
    c(
      utils::head(problems, max_n),
      sprintf("... and %d more", length(problems) - max_n)
    ),
    collapse = " | "
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
parse_markdown_table_row <- function(line) {
  line <- sub("^\\|", "", line)
  line <- sub("\\|[[:space:]]*$", "", line)
  fields <- trimws(strsplit(line, "\\|", fixed = FALSE)[[1]])
  if (length(fields) == 0L) {
    return(character())
  }
  fields
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
extract_markdown_link_target <- function(text) {
  match <- regexec("\\[[^]]+\\]\\(([^)]+)\\)", text, perl = TRUE)
  hit <- regmatches(text, match)[[1]]
  if (length(hit) < MARKDOWN_LINK_CAPTURE_LENGTH) {
    return(NA_character_)
  }
  hit[[MARKDOWN_LINK_CAPTURE_VALUE]]
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
frontmatter_scalar <- function(lines, key) {
  if (length(lines) < FRONTMATTER_MIN_LINES || !identical(lines[[1L]], "---")) {
    return(NA_character_)
  }
  end <- which(lines[-1L] == "---")
  if (length(end) == 0L) {
    return(NA_character_)
  }
  frontmatter <- lines[seq.int(2L, end[[1L]])]
  extract_yaml_scalar(frontmatter, key)
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
is_file_list_only_summary <- function(summary) {
  summary <- trimws(gsub("`", "", summary, fixed = TRUE))
  if (!nzchar(summary)) {
    return(FALSE)
  }
  summary <- trimws(sub("\\s*\\(\\+[0-9]+ more\\)\\s*$", "", summary))
  parts <- trimws(strsplit(summary, ",", fixed = TRUE)[[1]])
  if (length(parts) == 0L || any(!nzchar(parts))) {
    return(FALSE)
  }
  all(
    grepl("^[A-Za-z0-9_./ -]+$", parts) &
      grepl("\\.[A-Za-z0-9]+$", basename(parts))
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
make_tripwire_p27_sobj <- function() {
  cells <- paste0("cell", seq_len(12L))
  counts <- matrix(
    c(
      11,
      10,
      2,
      1,
      12,
      9,
      1,
      2,
      10,
      11,
      1,
      2,
      5,
      5,
      5,
      5,
      5,
      5,
      5,
      5,
      5,
      5,
      5,
      5
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("Cdkn1b", "Gapdh"), cells)
  )
  meta <- data.frame(
    # ANALYSIS_OK[sample-exclusion]: M1/M2/M3 are deliberately synthetic fixture identifiers; fixture cardinality and metadata fields are asserted below.
    Mouse = rep(c("M1", "M2", "M3"), each = 4L),
    Condition = rep(c("control", "estim", "control"), each = 4L),
    cluster = rep(c("1", "1", "2", "2"), times = 3L),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  # ANALYSIS_OK[warning-suppression]: Seurat fixture construction emits an expected assay-version warning; the resulting object is checked by downstream tripwire assertions.
  sobj <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    meta.data = meta
  ))
  sobj[["RNA"]] <- as(sobj[["RNA"]], Class = "Assay5")
  assay <- sobj[["RNA"]]
  SeuratObject::LayerData(assay, layer = "data") <- counts
  SeuratObject::LayerData(assay, layer = "pflog") <- counts
  sobj[["RNA"]] <- assay
  sobj
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
restore_random_seed <- function(seed_state, existed) {
  if (existed) {
    assign(".Random.seed", seed_state, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
analysis_table_annotation_dir <- function(root) {
  paths_env <- new.env(parent = .GlobalEnv)
  sys.source(file.path(root, "R", "paths.R"), envir = paths_env)
  file.path(get("TABLE_DIR", envir = paths_env, inherits = FALSE), "annotation")
}


MAX_PROBLEM_SUMMARY_ITEMS <- 8L
MARKDOWN_LINK_CAPTURE_LENGTH <- 2L
MARKDOWN_LINK_CAPTURE_VALUE <- 2L
FRONTMATTER_MIN_LINES <- 3L
STAGE_RECORD_WIDTH <- 3L
STAGE_RECORD_OFFSET <- 1L
MIN_EXPECTED_FINAL_STAGES <- 3L
MIN_MARKDOWN_TABLE_LINES <- 3L
MIN_METADATA_ROWS <- 2L
REGEX_CAPTURE_COUNT <- 2L
REGEX_CAPTURE_VALUE <- 2L
REGISTRY_HEADER_ROW_INDEX <- 1L:2L
PERMUTED_RUN_INDEX <- 2L

MAX_LINE_REFS <- 6L
YAML_LIST_CAPTURE_LENGTH <- 2L
YAML_LIST_CAPTURE_VALUE <- 2L


# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
line_refs <- function(lines, idx, max_n = MAX_LINE_REFS) {
  if (length(idx) == 0) {
    return("")
  }
  selected_idx <- utils::head(idx, max_n)
  paste(
    sprintf("L%d: %s", selected_idx, trimws(lines[selected_idx])),
    collapse = "; "
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
strip_inline_comment <- function(line) {
  # Good enough for the current scripts: ignore comment-only prose without trying
  # to parse quoted '#'. The tripwires below search code-like tokens.
  sub("^[[:space:]]*#.*$", "", line)
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
extract_yaml_list <- function(lines, key) {
  start <- grep(sprintf("^[[:space:]]*%s:[[:space:]]*$", key), lines)
  if (length(start) == 0) {
    return(character())
  }
  i <- start[[1]] + 1L
  out <- character()
  while (i <= length(lines)) {
    line <- lines[[i]]
    if (grepl("^[^[:space:]#][^:]*:", line)) {
      break
    }
    if (
      grepl("^[[:space:]]{0,2}[A-Za-z0-9_-]+:[[:space:]]*", line) &&
        !grepl("^[[:space:]]*-[[:space:]]*", line)
    ) {
      break
    }
    m <- regexec("^[[:space:]]*-[[:space:]]*\"?([^\"]+?)\"?[[:space:]]*$", line)
    hit <- regmatches(line, m)[[1]]
    if (length(hit) == YAML_LIST_CAPTURE_LENGTH) {
      out <- c(out, hit[[YAML_LIST_CAPTURE_VALUE]])
    }
    i <- i + 1L
  }
  out
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
extract_yaml_scalar <- function(lines, key) {
  idx <- grep(sprintf("^[[:space:]]*%s:[[:space:]]*", key), lines)
  if (length(idx) == 0) {
    return(NA_character_)
  }
  value <- sub(
    sprintf("^[[:space:]]*%s:[[:space:]]*", key),
    "",
    lines[[idx[[1]]]]
  )
  value <- sub("[[:space:]]+#.*$", "", value)
  value <- gsub("^\"|\"$", "", trimws(value))
  if (identical(value, "")) NA_character_ else value
}

# ANALYSIS_OK[file-freshness-tripwire]: mtime is the explicit freshness signal;
# report-values-freshness fails when rendered HTML predates source or figure targets.
mtime <- function(path) {
  info <- file.info(path)
  info$mtime[[1]]
}

resolve_figure_target <- function(path) {
  link_target <- Sys.readlink(path)
  if (!identical(link_target, "")) {
    if (!grepl("^/", link_target)) {
      link_target <- file.path(dirname(path), link_target)
    }
    return(normalizePath(link_target, winslash = "/", mustWork = FALSE))
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_branch_artifact_collision <- function(root) {
  # Scientific boundary: normalization and cell-cycle-filter branches must never
  # overwrite, share, or masquerade as the same clustering artifact. Persistent
  # Seurat artifact names must also avoid hyphens that Seurat rejects for
  # reduction names.
  slug <- "branch-artifact-collision"
  path <- file.path(root, "scripts", "04-cluster.R")
  if (!file.exists(path)) {
    return(fail(slug, "scripts/04-cluster.R is missing."))
  }
  lines <- read_text(path)
  code_lines <- strip_inline_comment(lines)
  code <- squash(code_lines)

  cc_start <- grep("\\bcc_tag[[:space:]]*<-", code_lines, perl = TRUE)
  cc_block <- integer()
  if (length(cc_start) > 0) {
    depth <- 0L
    for (i in seq.int(cc_start[[1]], length(code_lines))) {
      cc_block <- c(cc_block, i)
      chars <- strsplit(code_lines[[i]], "", fixed = TRUE)[[1]]
      depth <- depth +
        sum(chars %in% c("(", "{", "[")) -
        sum(chars %in% c(")", "}", "]"))
      if (i > cc_start[[1]] && depth <= 0L) {
        break
      }
    }
  }

  has_cc_tag <- length(cc_block) > 0 &&
    grepl("filtered_cell_cycle", squash(code_lines[cc_block]), fixed = TRUE)
  has_branch_tag <- grepl("\\bbranch_tag[[:space:]]*<-", code, perl = TRUE) &&
    grepl(
      "\\bbranch_tag[[:space:]]*<-.*\\bnorm\\b.*\\bcc_tag\\b|\\bbranch_tag[[:space:]]*<-.*\\bcc_tag\\b.*\\bnorm\\b",
      code,
      perl = TRUE
    )
  hyphenated_cc_tag <- cc_block[
    grepl(
      "\"(no-)?filter-cc\"|'(no-)?filter-cc'",
      code_lines[cc_block],
      perl = TRUE
    )
  ]
  has_branch_tag_guard <- grepl("\\bbranch_tag\\b", code, perl = TRUE) &&
    grepl("A-Za-z0-9_", code, fixed = TRUE) &&
    (grepl("^[A-Za-z0-9_]+$", code, fixed = TRUE) ||
      grepl("[^A-Za-z0-9_]", code, fixed = TRUE)) &&
    grepl("\\bgrepl[[:space:]]*\\(", code, perl = TRUE) &&
    grepl("\\bstop(ifnot)?[[:space:]]*\\(", code, perl = TRUE)

  persistent_name_line <- grepl(
    "sprintf\\([[:space:]]*\"(cluster_|umap_|%s_dims|.*\\.rds)",
    code_lines,
    perl = TRUE
  ) |
    grepl(
      "(prefix|out_tag|reduction_name|color_by)[[:space:]]*=",
      code_lines,
      perl = TRUE
    )
  norm_only <- which(
    persistent_name_line &
      grepl("\\bnorm\\b", code_lines, perl = TRUE) &
      !grepl("\\b(branch_tag|cc_tag)\\b", code_lines, perl = TRUE)
  )
  direct_branch_parts <- which(
    persistent_name_line &
      grepl("\\b(norm|cc_tag)\\b", code_lines, perl = TRUE) &
      !grepl("\\bbranch_tag\\b", code_lines, perl = TRUE)
  )
  literal_hyphen_persistent <- which(
    persistent_name_line &
      grepl("\"(no-)?filter-cc\"|'(no-)?filter-cc'", code_lines, perl = TRUE)
  )

  problems <- character()
  if (!has_cc_tag) {
    problems <- c(
      problems,
      "missing cc_tag derived from sobj@misc$preprocessing$filtered_cell_cycle"
    )
  }
  if (!has_branch_tag) {
    problems <- c(
      problems,
      "missing branch_tag that combines normalization and cell-cycle state"
    )
  }
  if (length(hyphenated_cc_tag) > 0) {
    problems <- c(
      problems,
      paste0(
        "cc_tag values feeding branch_tag are hyphenated; use Seurat-safe underscores: ",
        line_refs(lines, hyphenated_cc_tag)
      )
    )
  }
  if (!has_branch_tag_guard) {
    problems <- c(
      problems,
      "missing explicit branch_tag guard for allowed characters [A-Za-z0-9_]"
    )
  }
  if (length(norm_only) > 0) {
    problems <- c(
      problems,
      paste0(
        "persistent names still appear based on norm alone: ",
        line_refs(lines, norm_only)
      )
    )
  }
  if (length(direct_branch_parts) > 0) {
    problems <- c(
      problems,
      paste0(
        "persistent names should use the validated branch_tag instead of raw norm/cc_tag parts: ",
        line_refs(lines, direct_branch_parts)
      )
    )
  }
  if (length(literal_hyphen_persistent) > 0) {
    problems <- c(
      problems,
      paste0(
        "persistent Seurat names contain hyphenated cell-cycle tags: ",
        line_refs(lines, literal_hyphen_persistent)
      )
    )
  }

  if (length(problems) > 0) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }
  pass(
    slug,
    "cluster columns, UMAP reductions, clustree tags, and clustered RDS names use a validated normalization + cell-cycle branch_tag with Seurat-safe characters."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_cluster_wrapper_contract <- function(root) {
  # Operational boundary: the all-branch cluster wrapper must be an R
  # orchestrator with a non-executing preview path, not a shell loop depending
  # on exported state.
  slug <- "cluster-wrapper-contract"
  path <- file.path(root, "scripts", "04-cluster-all.R")
  if (!file.exists(path)) {
    return(fail(slug, "scripts/04-cluster-all.R is missing."))
  }
  lines <- read_text(path)
  code_lines <- strip_inline_comment(lines)
  code <- squash(code_lines)

  shell_loop <- grep(
    "\\bfor[[:space:]]+inp[[:space:]]+in\\b",
    code_lines,
    perl = TRUE
  )
  exported_current_dir <- grep(
    "(export[[:space:]]+CURRENT_OBJECT_DIR|\\$\\{?CURRENT_OBJECT_DIR\\}?)",
    code_lines,
    perl = TRUE
  )
  loads_package_constants <- grepl(
    "\\bdevtools::load_all[[:space:]]*\\(",
    code,
    perl = TRUE
  ) ||
    grepl(
      "\\blibrary[[:space:]]*\\([[:space:]]*ESPI[[:space:]]*\\)",
      code,
      perl = TRUE
    ) ||
    grepl(
      "\\brequireNamespace[[:space:]]*\\([[:space:]]*[\"']ESPI[\"']",
      code,
      perl = TRUE
    )
  uses_current_dir_in_r <- grepl(
    "\\bCURRENT_OBJECT_DIR\\b",
    code,
    perl = TRUE
  ) &&
    grepl("\\blist\\.files[[:space:]]*\\(", code, perl = TRUE) &&
    grepl("preprocess_.*\\.rds", code, perl = TRUE)
  has_preview_arg <- grepl(
    "--dry-run|dry[_\\.]run|preview",
    code,
    ignore.case = TRUE,
    perl = TRUE
  )
  prints_commands <- grepl(
    "\\b(writeLines|cat|message|print)[[:space:]]*\\(",
    code,
    perl = TRUE
  )
  builds_rscript_commands <- grepl("04-cluster.R", code, fixed = TRUE) &&
    grepl("\\b(Rscript|system2)[[:space:]]*\\(?", code, perl = TRUE)
  has_nonexecuting_preview <- has_preview_arg &&
    prints_commands &&
    builds_rscript_commands

  problems <- character()
  if (length(shell_loop) > 0) {
    problems <- c(
      problems,
      paste0(
        "cluster-all.R still contains a shell for-inp loop: ",
        line_refs(lines, shell_loop)
      )
    )
  }
  if (length(exported_current_dir) > 0) {
    problems <- c(
      problems,
      paste0(
        "cluster-all.R still references exported CURRENT_OBJECT_DIR shell text: ",
        line_refs(lines, exported_current_dir)
      )
    )
  }
  if (!loads_package_constants) {
    problems <- c(
      problems,
      "cluster-all.R does not load ESPI package constants in R"
    )
  }
  if (!uses_current_dir_in_r) {
    problems <- c(
      problems,
      "cluster-all.R does not discover preprocess_*.rds inputs from CURRENT_OBJECT_DIR in R"
    )
  }
  if (!has_nonexecuting_preview) {
    problems <- c(
      problems,
      "cluster-all.R lacks a dry-run/preview path that prints cluster-sobj.R commands without executing them"
    )
  }

  if (length(problems) > 0) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }
  pass(
    slug,
    "cluster-all.R loads package constants in R, discovers preprocess inputs without a shell loop, and exposes a non-executing command preview."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_cli_value_boundaries <- function(root) {
  # Operational boundary: value-bearing CLI flags must distinguish an absent
  # flag from a present flag with no value before any analysis artifacts are
  # regenerated.
  slug <- "cli-value-boundaries"
  rscript <- Sys.which("Rscript")
  if (identical(unname(rscript), "")) {
    return(fail(
      slug,
      "Rscript is unavailable, so CLI value-boundary behavior cannot be executed."
    ))
  }

  current_object_dir <- file.path(
    path.expand("~/Library/CloudStorage/Box-Box"),
    "megan_sc_data",
    "seurat_objects",
    "current"
  )
  preprocess_input <- file.path(
    current_object_dir,
    "preprocess_pflog_filter-cc.rds"
  )

  commands <- list(
    list(
      label = "04-cluster-all.R --dry-run --elbow-n",
      script = file.path(root, "scripts", "04-cluster-all.R"),
      args = c("--dry-run", "--elbow-n"),
      expected = "Missing value for --elbow-n"
    ),
    list(
      label = "03-preprocess.R --input --normalization pflog",
      script = file.path(root, "scripts", "03-preprocess.R"),
      args = c("--input", "--normalization", "pflog"),
      expected = "Missing value for --input"
    ),
    list(
      label = "03-preprocess.R invalid --input-source",
      script = file.path(root, "scripts", "03-preprocess.R"),
      args = c("--input-source", "invalid"),
      expected = "--input-source must be one of legacy or counts-qc"
    ),
    list(
      label = "03-preprocess-all.R mutually exclusive input options",
      script = file.path(root, "scripts", "03-preprocess-all.R"),
      args = c("--input", "object.rds", "--input-source", "legacy"),
      expected = "Use either --input or --input-source, not both."
    )
  )
  skipped <- character()
  if (file.exists(preprocess_input)) {
    commands <- c(
      commands,
      list(list(
        label = "04-cluster.R --input <valid-preprocess-object.rds> --elbow-n",
        script = file.path(root, "scripts", "04-cluster.R"),
        args = c("--input", preprocess_input, "--elbow-n"),
        expected = "Missing value for --elbow-n"
      ))
    )
  } else {
    skipped <- c(
      skipped,
      paste(
        "04-cluster.R valueless --elbow-n subcase skipped because",
        preprocess_input,
        "does not exist"
      )
    )
  }

  failures <- character()
  for (command in commands) {
    if (!file.exists(command$script)) {
      failures <- c(failures, paste(command$script, "is missing"))
      next
    }
    # ANALYSIS_OK[warning-suppression]: nonzero exits are the expected signal
    # in this malformed-CLI tripwire; captured output is checked below.
    output <- tryCatch(
      suppressWarnings(system2(
        rscript,
        c(command$script, command$args),
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(e) structure(conditionMessage(e), status = 1L)
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }
    text <- paste(output, collapse = "\n")

    if (identical(as.integer(status), 0L)) {
      failures <- c(
        failures,
        paste("unexpected success for", command$label)
      )
    } else if (!grepl(command$expected, text, fixed = TRUE)) {
      failures <- c(
        failures,
        paste(
          "missing expected error for",
          command$label,
          ":",
          command$expected
        )
      )
    }
  }

  if (length(failures) > 0) {
    message <- paste(failures, collapse = " | ")
    if (length(skipped) > 0) {
      message <- paste(message, paste(skipped, collapse = " | "))
    }
    return(fail(slug, message))
  }

  message <- "Value-bearing CLI flags fail when present without values."
  if (length(skipped) > 0) {
    message <- paste(message, paste(skipped, collapse = " | "))
  }
  pass(slug, message)
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_pipeline_dry_run_contract <- function(root) {
  # Public boundary: dry-run exposes a complete, safe, line-oriented plan.
  slug <- "pipeline-dry-run-contract"
  rscript <- Sys.which("Rscript")
  if (identical(unname(rscript), "")) {
    return(fail(
      slug,
      "Rscript is unavailable, so the public pipeline dry-run cannot be executed."
    ))
  }

  script <- file.path(root, "scripts", "run-pipeline.R")
  if (!file.exists(script)) {
    return(fail(slug, "scripts/run-pipeline.R is missing."))
  }

  counts_qc_stages <- c(
    "process-counts",
    "qc-filtering",
    "preprocess-source",
    "cluster-source-log1p-no-filter-cc",
    "cluster-source-log1p-filter-cc",
    "cluster-source-pflog-no-filter-cc",
    "cluster-source-pflog-filter-cc",
    "summarize-source",
    "select-mg",
    "cluster-mg-no-filter-cc",
    "cluster-mg-filter-cc",
    "summarize-mg",
    "marker-heatmap-source",
    "marker-heatmap-mg-no-filter-cc",
    "marker-heatmap-mg-filter-cc",
    "module-heatmap-source",
    "module-heatmap-mg-no-filter-cc",
    "module-heatmap-mg-filter-cc",
    "mg-figures-no-filter-cc",
    "mg-figures-filter-cc",
    "mg-markers",
    "mg-de",
    "render-notebook",
    "tripwires"
  )
  existing_source_stages <- counts_qc_stages[-c(1L, 2L)]
  explicit_input <- file.path(root, "data", "explicit-input.seurat.rds")

  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  make_header <- function(input_source, overwrite, stages) {
    c(
      "mode: dry-run",
      paste0("input_source: ", input_source),
      paste0("overwrite: ", overwrite),
      "source_cluster_column: cluster_pflog_no_filter_cc_dims30_res0.3",
      "mg_cluster_column: cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5",
      "mg_pca_dims: 50",
      paste0("first_stage: ", stages[[1L]]),
      "final_stage: tripwires",
      paste0("stage_count: ", length(stages))
    )
  }
  successful_invocations <- list(
    list(
      label = "counts-qc",
      args = c("--dry-run"),
      input_source = "counts-qc",
      overwrite = FALSE,
      stages = counts_qc_stages
    ),
    list(
      label = "legacy",
      args = c("--dry-run", "--input-source", "legacy"),
      input_source = "legacy",
      overwrite = FALSE,
      stages = existing_source_stages
    ),
    list(
      label = "explicit",
      args = c("--dry-run", "--input", shQuote(explicit_input)),
      input_source = "explicit",
      overwrite = FALSE,
      stages = existing_source_stages
    ),
    list(
      label = "counts-qc overwrite",
      args = c("--dry-run", "--overwrite"),
      input_source = "counts-qc",
      overwrite = TRUE,
      stages = counts_qc_stages
    )
  )
  failed_invocations <- list(
    list(
      label = "missing input source value (--dry-run --input-source)",
      args = c("--dry-run", "--input-source"),
      expected = "Missing value for --input-source."
    ),
    list(
      label = "missing input value (--dry-run --input)",
      args = c("--dry-run", "--input"),
      expected = "Missing value for --input."
    ),
    list(
      label = "invalid input source (--dry-run --input-source invalid)",
      args = c("--dry-run", "--input-source", "invalid"),
      expected = "--input-source must be one of counts-qc or legacy."
    ),
    list(
      label = "invalid input source (--input-source invalid)",
      args = c("--input-source", "invalid"),
      expected = "--input-source must be one of counts-qc or legacy."
    ),
    list(
      label = "mutually exclusive input options",
      args = c(
        "--dry-run",
        "--input",
        shQuote(explicit_input),
        "--input-source",
        "legacy"
      ),
      expected = "Use either --input or --input-source, not both."
    ),
    list(
      label = "unknown flag (--dry-run --unknown)",
      args = c("--dry-run", "--unknown"),
      expected = "Unknown argument: --unknown."
    )
  )
  run_pipeline <- function(args) {
    # ANALYSIS_OK[warning-suppression]: system2() emits expected command-not-found diagnostics while probing CLI failures; captured output and status are asserted below.
    output <- tryCatch(
      suppressWarnings(system2(
        rscript,
        c(shQuote(script), args),
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(e) structure(conditionMessage(e), status = 1L)
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }
    list(status = as.integer(status), output = unname(as.character(output)))
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  format_output <- function(output) {
    paste(shQuote(output), collapse = ", ")
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  output_at <- function(output, position) {
    if (position > length(output)) NA_character_ else output[[position]]
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  command_has <- function(command, value) {
    !is.na(command) && grepl(value, command, fixed = TRUE)
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  path_has_suffix <- function(path, suffix) {
    if (is.na(path)) {
      return(FALSE)
    }
    path <- gsub("\\\\", "/", path)
    suffix <- gsub("\\\\", "/", suffix)
    identical(path, suffix) || endsWith(path, paste0("/", suffix))
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  stage_problem <- function(variant, stage, expected, actual) {
    paste0(
      variant,
      ": stage ",
      stage,
      ": expected ",
      expected,
      "; got ",
      shQuote(actual)
    )
  }

  problems <- character()
  for (invocation in successful_invocations) {
    result <- run_pipeline(invocation$args)
    if (!identical(result$status, 0L)) {
      problems <- c(
        problems,
        sprintf(
          "%s: expected exit status 0, got %d",
          invocation$label,
          result$status
        )
      )
    }

    header <- make_header(
      invocation$input_source,
      tolower(as.character(invocation$overwrite)),
      invocation$stages
    )
    if (
      length(result$output) < length(header) ||
        !identical(result$output[seq_along(header)], header)
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": header mismatch; expected ",
          format_output(header),
          "; got ",
          format_output(utils::head(result$output, length(header)))
        )
      )
    }

    expected_line_count <- length(header) + 3L * length(invocation$stages)
    if (!identical(length(result$output), expected_line_count)) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": expected ",
          expected_line_count,
          " contract lines for ",
          length(invocation$stages),
          " stages; got ",
          length(result$output)
        )
      )
    }

    stage_records <- vector("list", length(invocation$stages))
    names(stage_records) <- invocation$stages
    observed_stages <- character(length(invocation$stages))
    for (stage_index in seq_along(invocation$stages)) {
      stage <- invocation$stages[[stage_index]]
      record_start <- length(header) +
        (stage_index - STAGE_RECORD_OFFSET) * STAGE_RECORD_WIDTH +
        STAGE_RECORD_OFFSET
      actual_stage <- output_at(result$output, record_start)
      observed_stages[[stage_index]] <- sub("^stage: ", "", actual_stage)
      expected_stage <- paste0("stage: ", stage)
      if (!identical(actual_stage, expected_stage)) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            stage,
            paste0("stage line ", shQuote(expected_stage)),
            actual_stage
          )
        )
      }

      command <- output_at(result$output, record_start + 1L)
      if (
        is.na(command) ||
          !grepl("^command: ('[^']*')( ('[^']*'))*$", command, perl = TRUE)
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            stage,
            "one shell-quoted command line",
            command
          )
        )
      }

      expects <- output_at(result$output, record_start + 2L)
      if (
        is.na(expects) ||
          !grepl("^expects: [^;[:space:]]+(;[^;[:space:]]+)*$", expects)
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            stage,
            "one-or-more semicolon-separated expected output paths",
            expects
          )
        )
      }
      stage_records[[stage]] <- list(command = command, expects = expects)
    }
    # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
    stage_output_paths <- function(record) {
      if (is.null(record) || is.na(record$expects)) {
        return(character())
      }
      strsplit(sub("^expects: ", "", record$expects), ";", fixed = TRUE)[[1L]]
    }
    expected_count_output <- "data/input/sobj_raw.rds"
    expected_qc_outputs <- c(
      "seurat_objects/input/sobj_raw_with_qc.rds",
      "seurat_objects/input/sobj_qc_filtered.rds"
    )
    if (identical(invocation$input_source, "counts-qc")) {
      counts_outputs <- stage_output_paths(stage_records[["process-counts"]])
      if (
        !identical(length(counts_outputs), 1L) ||
          !endsWith(counts_outputs, expected_count_output)
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            "process-counts",
            paste0(
              "exactly one DATA_ROOT_DIR/",
              expected_count_output,
              " output"
            ),
            stage_records[["process-counts"]]$expects
          )
        )
      }

      qc_outputs <- stage_output_paths(stage_records[["qc-filtering"]])
      if (
        !identical(length(qc_outputs), length(expected_qc_outputs)) ||
          !all(vapply(
            expected_qc_outputs,
            function(output) any(endsWith(qc_outputs, output)),
            logical(1)
          ))
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            "qc-filtering",
            paste0(
              "both INPUT_OBJECT_DIR outputs ",
              paste(expected_qc_outputs, collapse = ", ")
            ),
            stage_records[["qc-filtering"]]$expects
          )
        )
      }
    }

    preprocess_command <- stage_records[["preprocess-source"]]$command
    expected_source_option <- if (
      identical(invocation$input_source, "explicit")
    ) {
      paste0("'--input' ", shQuote(explicit_input))
    } else {
      paste0("'--input-source' '", invocation$input_source, "'")
    }
    if (
      !command_has(preprocess_command, "'scripts/03-preprocess-all.R'") ||
        !command_has(preprocess_command, expected_source_option)
    ) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "preprocess-source",
          paste0(
            "scripts/03-preprocess-all.R command containing ",
            shQuote(expected_source_option)
          ),
          preprocess_command
        )
      )
    }

    source_clusters <- c(
      "cluster-source-log1p-no-filter-cc" = "preprocess_log1p_no-filter-cc.rds",
      "cluster-source-log1p-filter-cc" = "preprocess_log1p_filter-cc.rds",
      "cluster-source-pflog-no-filter-cc" = "preprocess_pflog_no-filter-cc.rds",
      "cluster-source-pflog-filter-cc" = "preprocess_pflog_filter-cc.rds"
    )
    source_cluster_outputs <- c(
      "cluster-source-log1p-no-filter-cc" = "cluster_log1p_no_filter_cc_elbow20.rds",
      "cluster-source-log1p-filter-cc" = "cluster_log1p_filter_cc_elbow20.rds",
      "cluster-source-pflog-no-filter-cc" = "cluster_pflog_no_filter_cc_elbow20.rds",
      "cluster-source-pflog-filter-cc" = "cluster_pflog_filter_cc_elbow20.rds"
    )
    current_object_suffix <- file.path("seurat_objects", "current")
    preprocess_expects <- stage_records[["preprocess-source"]]$expects
    preprocess_outputs <- if (is.na(preprocess_expects)) {
      character()
    } else {
      strsplit(sub("^expects: ", "", preprocess_expects), ";", fixed = TRUE)[[
        1L
      ]]
    }
    preprocess_matches <- length(preprocess_outputs) ==
      length(source_clusters) &&
      all(vapply(
        seq_along(source_clusters),
        function(index) {
          path_has_suffix(
            preprocess_outputs[[index]],
            file.path(current_object_suffix, unname(source_clusters)[[index]])
          )
        },
        logical(1)
      ))
    if (!isTRUE(preprocess_matches)) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "preprocess-source",
          paste(
            "exactly the four source preprocess outputs under",
            current_object_suffix,
            paste(unname(source_clusters), collapse = ", ")
          ),
          preprocess_expects
        )
      )
    }
    source_summary_outputs <- file.path(
      c(
        "tables",
        "tables",
        "tables",
        "figures",
        "figures",
        "figures",
        "figures"
      ),
      c(
        "cluster/cluster_grid_summary.tsv",
        "cluster/cluster_grid_stability_summary.tsv",
        "cluster/cluster_grid_pairwise_stability.tsv",
        "cluster/cluster_grid_clustree_12_panel.png",
        "cluster/cluster_grid_clustree_12_panel.pdf",
        "cluster/umap_resolution_sweep_pflog_filter_cc_dims50.png",
        "cluster/umap_resolution_sweep_pflog_filter_cc_dims50.pdf"
      )
    )
    cluster_execution_options <- c(
      "'--elbow-n' '20'",
      "'--extra-dims' '30,50'",
      "'--resolutions' '0.3,0.5,0.8'"
    )
    cluster_candidate_dims <- c(20L, 30L, 50L)
    cluster_candidate_resolution_tags <- gsub(
      "[^A-Za-z0-9_-]",
      "_",
      c("0.3", "0.5", "0.8")
    )
    # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
    cluster_notebook_umap_outputs <- function(branch) {
      unlist(
        lapply(
          cluster_candidate_dims,
          function(dims) {
            file.path(
              "notebook",
              "figures",
              sprintf(
                "umap_%s_dims%d_by_cluster_%s_dims%d_res%s.png",
                branch,
                dims,
                branch,
                dims,
                cluster_candidate_resolution_tags
              )
            )
          }
        ),
        use.names = FALSE
      )
    }
    # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
    cluster_stage_expected_outputs <- function(rds_path) {
      branch <- sub(
        "^cluster_",
        "",
        sub("_elbow20\\.rds$", "", basename(rds_path))
      )
      c(
        rds_path,
        cluster_notebook_umap_outputs(branch)
      )
    }
    # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
    cluster_stage_expected_description <- function(rds_path, branch_label) {
      paste(
        "exactly one clustered",
        branch_label,
        "RDS plus nine notebook UMAP PNG links",
        paste(cluster_stage_expected_outputs(rds_path), collapse = ", ")
      )
    }
    # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
    stage_outputs_match <- function(record, expected_suffixes) {
      actual_outputs <- stage_output_paths(record)
      length(actual_outputs) == length(expected_suffixes) &&
        all(vapply(
          expected_suffixes,
          function(suffix) {
            sum(vapply(
              actual_outputs,
              path_has_suffix,
              logical(1),
              suffix = suffix
            )) ==
              1L
          },
          logical(1)
        ))
    }
    for (stage in names(source_clusters)) {
      command <- stage_records[[stage]]$command
      expected_input <- source_clusters[[stage]]
      if (
        !command_has(command, "'scripts/04-cluster.R'") ||
          !command_has(command, expected_input)
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            stage,
            paste0(
              "scripts/04-cluster.R command for ",
              shQuote(expected_input)
            ),
            command
          )
        )
      }
      for (expected_option in cluster_execution_options) {
        if (!command_has(command, expected_option)) {
          problems <- c(
            problems,
            stage_problem(
              invocation$label,
              stage,
              paste0("command containing ", shQuote(expected_option)),
              command
            )
          )
        }
      }
      if (
        !stage_outputs_match(
          stage_records[[stage]],
          cluster_stage_expected_outputs(
            file.path(current_object_suffix, source_cluster_outputs[[stage]])
          )
        )
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            stage,
            cluster_stage_expected_description(
              file.path(current_object_suffix, source_cluster_outputs[[stage]]),
              "source"
            ),
            stage_records[[stage]]$expects
          )
        )
      }
    }

    summary_command <- stage_records[["summarize-source"]]$command
    if (!command_has(summary_command, "'scripts/05-summarize-clusters.R'")) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "summarize-source",
          "scripts/05-summarize-clusters.R command",
          summary_command
        )
      )
    }
    summary_expects <- stage_output_paths(stage_records[["summarize-source"]])
    summary_matches <- length(summary_expects) ==
      length(source_summary_outputs) &&
      all(vapply(
        source_summary_outputs,
        function(expected) {
          any(vapply(
            summary_expects,
            path_has_suffix,
            logical(1),
            suffix = expected
          ))
        },
        logical(1)
      ))
    if (!isTRUE(summary_matches)) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "summarize-source",
          paste(
            "all source summary artifacts",
            paste(source_summary_outputs, collapse = ", ")
          ),
          stage_records[["summarize-source"]]$expects
        )
      )
    }

    select_mg_command <- stage_records[["select-mg"]]$command
    select_mg_options <- c(
      "'--cluster-column' 'cluster_pflog_no_filter_cc_dims30_res0.3'",
      "'--dims' '50'"
    )
    for (expected_option in select_mg_options) {
      if (!command_has(select_mg_command, expected_option)) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            "select-mg",
            paste0("command containing ", shQuote(expected_option)),
            select_mg_command
          )
        )
      }
    }

    select_mg_outputs <- c(
      "seurat_objects/current/preprocess_pflog_mg_selected_no-filter-cc.rds",
      "seurat_objects/current/preprocess_pflog_mg_selected_filter-cc.rds",
      "tables/mg_selected/mg_selected_cluster_selection.tsv",
      "figures/mg_selected/mg_selected_cluster_selection_diagnostics.png",
      "figures/mg_selected/mg_selected_cluster_selection_diagnostics.pdf",
      "figures/mg_selected/elbow_pflog_mg_selected_no_filter_cc.png",
      "figures/mg_selected/elbow_pflog_mg_selected_no_filter_cc.pdf",
      "figures/mg_selected/elbow_pflog_mg_selected_filter_cc.png",
      "figures/mg_selected/elbow_pflog_mg_selected_filter_cc.pdf"
    )
    if (
      !stage_outputs_match(
        stage_records[["select-mg"]],
        select_mg_outputs
      )
    ) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "select-mg",
          paste(
            "exactly the selection RDS, table, diagnostics, and elbow artifacts",
            paste(select_mg_outputs, collapse = ", ")
          ),
          stage_records[["select-mg"]]$expects
        )
      )
    }
    figure_stage_contracts <- list(
      "marker-heatmap-source" = list(
        script = "scripts/06-plot-marker-heatmap.R",
        options = c("'--dims' '30'", "'--resolution' '0.3'"),
        outputs = c(
          file.path(
            "figures/annotation",
            paste0(
              "cell_type_marker_heatmap_pflog_pflog_no_filter_cc_cells_dims30_res0.3",
              c(".png", ".pdf")
            )
          ),
          file.path(
            "notebook/figures",
            "cell_type_marker_heatmap_pflog_pflog_no_filter_cc_cells_dims30_res0.3.png"
          )
        )
      ),
      "marker-heatmap-mg-no-filter-cc" = list(
        script = "scripts/06-plot-marker-heatmap.R",
        options = c("'--dims' '20'", "'--resolution' '0.5'"),
        outputs = c(
          file.path(
            "figures/annotation",
            paste0(
              "cell_type_marker_heatmap_pflog_pflog_mg_selected_no_filter_cc_cells_dims20_res0.5",
              c(".png", ".pdf")
            )
          ),
          file.path(
            "notebook/figures",
            "cell_type_marker_heatmap_pflog_pflog_mg_selected_no_filter_cc_cells_dims20_res0.5.png"
          )
        )
      ),
      "marker-heatmap-mg-filter-cc" = list(
        script = "scripts/06-plot-marker-heatmap.R",
        options = c("'--dims' '20'", "'--resolution' '0.5'"),
        outputs = c(
          file.path(
            "figures/annotation",
            paste0(
              "cell_type_marker_heatmap_pflog_pflog_mg_selected_filter_cc_cells_dims20_res0.5",
              c(".png", ".pdf")
            )
          ),
          file.path(
            "notebook/figures",
            "cell_type_marker_heatmap_pflog_pflog_mg_selected_filter_cc_cells_dims20_res0.5.png"
          )
        )
      ),
      "module-heatmap-source" = list(
        script = "scripts/10-plot-cluster-marker-heatmaps.R",
        options = c("'--dims' '30'", "'--resolution' '0.3'"),
        outputs = c(
          file.path(
            "figures/annotation",
            paste0(
              "cell_type_module_p27_heatmap_pflog_pflog_no_filter_cc_dims30_res0.3",
              c(".png", ".pdf")
            )
          ),
          file.path(
            "tables/annotation",
            paste0(
              "cell_type_module_p27_heatmap_pflog_pflog_no_filter_cc_dims30_res0.3",
              c("_module_scores.tsv", "_p27_enrichment.tsv")
            )
          ),
          file.path(
            "notebook/figures",
            "cell_type_module_p27_heatmap_pflog_pflog_no_filter_cc_dims30_res0.3.png"
          )
        )
      ),
      "module-heatmap-mg-no-filter-cc" = list(
        script = "scripts/10-plot-cluster-marker-heatmaps.R",
        options = c("'--dims' '20'", "'--resolution' '0.5'"),
        outputs = c(
          file.path(
            "figures/annotation",
            paste0(
              "cell_type_module_p27_heatmap_pflog_pflog_mg_selected_no_filter_cc_dims20_res0.5",
              c(".png", ".pdf")
            )
          ),
          file.path(
            "tables/annotation",
            paste0(
              "cell_type_module_p27_heatmap_pflog_pflog_mg_selected_no_filter_cc_dims20_res0.5",
              c("_module_scores.tsv", "_p27_enrichment.tsv")
            )
          ),
          file.path(
            "notebook/figures",
            "cell_type_module_p27_heatmap_pflog_pflog_mg_selected_no_filter_cc_dims20_res0.5.png"
          )
        )
      ),
      "module-heatmap-mg-filter-cc" = list(
        script = "scripts/10-plot-cluster-marker-heatmaps.R",
        options = c("'--dims' '20'", "'--resolution' '0.5'"),
        outputs = c(
          file.path(
            "figures/annotation",
            paste0(
              "cell_type_module_p27_heatmap_pflog_pflog_mg_selected_filter_cc_dims20_res0.5",
              c(".png", ".pdf")
            )
          ),
          file.path(
            "tables/annotation",
            paste0(
              "cell_type_module_p27_heatmap_pflog_pflog_mg_selected_filter_cc_dims20_res0.5",
              c("_module_scores.tsv", "_p27_enrichment.tsv")
            )
          ),
          file.path(
            "notebook/figures",
            "cell_type_module_p27_heatmap_pflog_pflog_mg_selected_filter_cc_dims20_res0.5.png"
          )
        )
      ),
      "mg-figures-no-filter-cc" = list(
        script = "scripts/09-plot-mg-figures.R",
        options = c(
          "'--elbow-n' '20'",
          "'--dims' '20'",
          "'--resolution' '0.5'"
        ),
        outputs = c(
          file.path(
            "figures/mg_selected",
            paste0(
              c(
                "mg_selected_cluster_umap_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_condition_umap_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_feature_umap_pflog_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_ascl1_hes6_coexpression_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_cluster_abundance_enrichment_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_cluster_proportion_by_mouse_pflog_mg_selected_no_filter_cc_dims20_res0.5"
              ),
              rep(c(".png", ".pdf"), each = 6L)
            )
          ),
          file.path(
            "tables/mg_selected",
            paste0(
              c(
                "mg_selected_cluster_abundance_enrichment_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_cluster_proportion_randomization_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_sample_cluster_proportions_pflog_mg_selected_no_filter_cc_dims20_res0.5"
              ),
              ".tsv"
            )
          ),
          file.path(
            "notebook/figures",
            paste0(
              c(
                "mg_selected_cluster_umap_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_condition_umap_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_feature_umap_pflog_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_ascl1_hes6_coexpression_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_cluster_abundance_enrichment_pflog_mg_selected_no_filter_cc_dims20_res0.5",
                "mg_selected_cluster_proportion_by_mouse_pflog_mg_selected_no_filter_cc_dims20_res0.5"
              ),
              ".png"
            )
          )
        )
      ),
      "mg-figures-filter-cc" = list(
        script = "scripts/09-plot-mg-figures.R",
        options = c(
          "'--elbow-n' '20'",
          "'--dims' '20'",
          "'--resolution' '0.5'"
        ),
        outputs = c(
          file.path(
            "figures/mg_selected",
            paste0(
              c(
                "mg_selected_cluster_umap_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_condition_umap_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_feature_umap_pflog_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_ascl1_hes6_coexpression_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_cluster_abundance_enrichment_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_cluster_proportion_by_mouse_pflog_mg_selected_filter_cc_dims20_res0.5"
              ),
              rep(c(".png", ".pdf"), each = 6L)
            )
          ),
          file.path(
            "tables/mg_selected",
            paste0(
              c(
                "mg_selected_cluster_abundance_enrichment_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_cluster_proportion_randomization_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_sample_cluster_proportions_pflog_mg_selected_filter_cc_dims20_res0.5"
              ),
              ".tsv"
            )
          ),
          file.path(
            "notebook/figures",
            paste0(
              c(
                "mg_selected_cluster_umap_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_condition_umap_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_feature_umap_pflog_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_ascl1_hes6_coexpression_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_cluster_abundance_enrichment_pflog_mg_selected_filter_cc_dims20_res0.5",
                "mg_selected_cluster_proportion_by_mouse_pflog_mg_selected_filter_cc_dims20_res0.5"
              ),
              ".png"
            )
          )
        )
      )
    )
    for (stage in names(figure_stage_contracts)) {
      contract <- figure_stage_contracts[[stage]]
      command <- stage_records[[stage]]$command
      for (expected_token in c(
        paste0("'", contract$script, "'"),
        contract$options
      )) {
        if (!command_has(command, expected_token)) {
          problems <- c(
            problems,
            stage_problem(
              invocation$label,
              stage,
              paste0("command containing ", shQuote(expected_token)),
              command
            )
          )
        }
      }
      if (
        !stage_outputs_match(
          stage_records[[stage]],
          contract$outputs
        )
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            stage,
            paste(
              "exactly the chosen-dimension/resolution artifacts",
              paste(contract$outputs, collapse = ", ")
            ),
            stage_records[[stage]]$expects
          )
        )
      }
    }

    mg_cluster_outputs <- c(
      "cluster-mg-no-filter-cc" = "seurat_objects/current/cluster_pflog_mg_selected_no_filter_cc_elbow20.rds",
      "cluster-mg-filter-cc" = "seurat_objects/current/cluster_pflog_mg_selected_filter_cc_elbow20.rds"
    )
    for (stage in names(mg_cluster_outputs)) {
      command <- stage_records[[stage]]$command
      for (expected_option in cluster_execution_options) {
        if (!command_has(command, expected_option)) {
          problems <- c(
            problems,
            stage_problem(
              invocation$label,
              stage,
              paste0("command containing ", shQuote(expected_option)),
              command
            )
          )
        }
      }
      if (
        !stage_outputs_match(
          stage_records[[stage]],
          cluster_stage_expected_outputs(mg_cluster_outputs[[stage]])
        )
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            stage,
            cluster_stage_expected_description(
              mg_cluster_outputs[[stage]],
              "MG"
            ),
            stage_records[[stage]]$expects
          )
        )
      }
    }

    summarize_mg_command <- stage_records[["summarize-mg"]]$command
    if (
      !command_has(
        summarize_mg_command,
        "'scripts/08-summarize-mg-clusters.R'"
      ) ||
        !command_has(summarize_mg_command, "'--elbow-n' '20'")
    ) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "summarize-mg",
          "scripts/08-summarize-mg-clusters.R command containing '--elbow-n' '20'",
          summarize_mg_command
        )
      )
    }
    summarize_mg_outputs <- c(
      "tables/mg_selected/mg_selected_cluster_grid_summary.tsv",
      "figures/mg_selected/mg_selected_umap_resolution_sweep_pflog_mg_selected_no_filter_cc_dims50.png",
      "figures/mg_selected/mg_selected_umap_resolution_sweep_pflog_mg_selected_no_filter_cc_dims50.pdf",
      "figures/mg_selected/mg_selected_umap_resolution_sweep_pflog_mg_selected_filter_cc_dims50.png",
      "figures/mg_selected/mg_selected_umap_resolution_sweep_pflog_mg_selected_filter_cc_dims50.pdf"
    )
    if (
      !stage_outputs_match(
        stage_records[["summarize-mg"]],
        summarize_mg_outputs
      )
    ) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "summarize-mg",
          paste(
            "exactly the MG grid summary and no-filter/filter-CC dims50 sweeps",
            paste(summarize_mg_outputs, collapse = ", ")
          ),
          stage_records[["summarize-mg"]]$expects
        )
      )
    }

    expected_marker_outputs <- c(
      "tables/mg_selected/find_all_markers_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv",
      "tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv",
      "tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv",
      "tables/mg_selected/find_all_markers_identity_map_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv",
      "figures/mg_selected/mg_selected_cluster_marker_dotplot_data_pflog_mg_selected_no_filter_cc_dims20_res0.5_top5.png",
      "figures/mg_selected/mg_selected_cluster_marker_dotplot_data_pflog_mg_selected_no_filter_cc_dims20_res0.5_top5.pdf",
      "notebook/figures/mg_selected_cluster_marker_dotplot_data_pflog_mg_selected_no_filter_cc_dims20_res0.5_top5.png"
    )
    expected_de_outputs <- c(
      "degs/mg_selected/pseudobulk_sample_summary.tsv",
      "degs/mg_selected/design_summary.tsv",
      "degs/mg_selected/deseq2_full_results.tsv",
      "degs/mg_selected/deseq2_significant_degs.tsv",
      "degs/mg_selected/deseq2_marker_overlap.tsv",
      "degs/mg_selected/deseq2_paired_sensitivity_full_results.tsv",
      "degs/mg_selected/deseq2_paired_sensitivity_significant_degs.tsv",
      "degs/mg_selected/deseq2_paired_sensitivity_marker_overlap.tsv",
      "degs/mg_selected/numbers.json",
      "enrichment/mg_selected/go_bp_ora_up.tsv",
      "enrichment/mg_selected/go_bp_ora_down.tsv",
      "enrichment/mg_selected/go_bp_gsea.tsv",
      "enrichment/mg_selected/go_bp_gsea_symbol_entrez_mapping.tsv",
      "enrichment/mg_selected/go_bp_ora_up_simplified.tsv",
      "enrichment/mg_selected/go_bp_ora_down_simplified.tsv",
      "enrichment/mg_selected/go_bp_gsea_simplified.tsv",
      "enrichment/mg_selected/go_bp_ora_up_bayes_simplified.tsv",
      "enrichment/mg_selected/go_bp_ora_down_bayes_simplified.tsv",
      "figures/mg_selected/mg_selected_go_ora_up_dotplot.png",
      "figures/mg_selected/mg_selected_go_ora_up_dotplot.pdf",
      "figures/mg_selected/mg_selected_go_ora_down_dotplot.png",
      "figures/mg_selected/mg_selected_go_ora_down_dotplot.pdf",
      "figures/mg_selected/mg_selected_go_gsea_dotplot.png",
      "figures/mg_selected/mg_selected_go_gsea_dotplot.pdf",
      "figures/mg_selected/mg_selected_go_ora_up_bayes_dotplot.png",
      "figures/mg_selected/mg_selected_go_ora_up_bayes_dotplot.pdf",
      "figures/mg_selected/mg_selected_go_ora_down_bayes_dotplot.png",
      "figures/mg_selected/mg_selected_go_ora_down_bayes_dotplot.pdf",
      "notebook/figures/mg_selected_go_ora_up_dotplot.png",
      "notebook/figures/mg_selected_go_ora_down_dotplot.png",
      "notebook/figures/mg_selected_go_gsea_dotplot.png",
      "notebook/figures/mg_selected_go_ora_up_bayes_dotplot.png",
      "notebook/figures/mg_selected_go_ora_down_bayes_dotplot.png",
      "figures/mg_selected/mg_selected_de_volcano.png",
      "figures/mg_selected/mg_selected_de_volcano.pdf",
      "notebook/figures/mg_selected_de_volcano.png"
    )

    mg_marker_command <- stage_records[["mg-markers"]]$command
    expected_marker_options <- c(
      "'scripts/11-find-mg-markers.R'",
      "cluster_pflog_mg_selected_no_filter_cc_elbow20.rds",
      "'--branch-tag' 'pflog_mg_selected_no_filter_cc'",
      "'--dims' '20'",
      "'--resolution' '0.5'",
      "'--layer' 'data'",
      "'--counts-layer' 'counts'",
      "'--confirm-no-merge'"
    )
    for (expected_option in expected_marker_options) {
      if (!command_has(mg_marker_command, expected_option)) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            "mg-markers",
            paste0("command containing ", shQuote(expected_option)),
            mg_marker_command
          )
        )
      }
    }
    if (
      !stage_outputs_match(
        stage_records[["mg-markers"]],
        expected_marker_outputs
      )
    ) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "mg-markers",
          paste(
            "all seven protected MG marker outputs",
            paste(expected_marker_outputs, collapse = ", ")
          ),
          stage_records[["mg-markers"]]$expects
        )
      )
    }

    mg_de_command <- stage_records[["mg-de"]]$command
    expected_de_options <- c(
      "'scripts/12-run-mg-de.R'",
      "cluster_pflog_mg_selected_no_filter_cc_elbow20.rds",
      "'--cluster-column' 'cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5'",
      "'--counts-layer' 'counts'",
      "'--lfc-shrink-type' 'apeglm'"
    )
    for (expected_option in expected_de_options) {
      if (!command_has(mg_de_command, expected_option)) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            "mg-de",
            paste0("command containing ", shQuote(expected_option)),
            mg_de_command
          )
        )
      }
    }
    if (
      !stage_outputs_match(
        stage_records[["mg-de"]],
        expected_de_outputs
      )
    ) {
      problems <- c(
        problems,
        stage_problem(
          invocation$label,
          "mg-de",
          paste(
            "all protected MG DE outputs",
            paste(expected_de_outputs, collapse = ", ")
          ),
          stage_records[["mg-de"]]$expects
        )
      )
    }

    overwrite_stages <- names(stage_records)[vapply(
      stage_records,
      function(record) command_has(record$command, "'--overwrite'"),
      logical(1)
    )]
    expected_overwrite_stages <- if (isTRUE(invocation$overwrite)) {
      c("mg-markers", "mg-de")
    } else {
      character()
    }
    if (!identical(overwrite_stages, expected_overwrite_stages)) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": --overwrite stage mismatch; expected ",
          paste(expected_overwrite_stages, collapse = ", "),
          "; got ",
          paste(overwrite_stages, collapse = ", ")
        )
      )
    }

    marker_index <- match("mg-markers", observed_stages)
    de_index <- match("mg-de", observed_stages)
    figure_indices <- match(names(figure_stage_contracts), observed_stages)
    last_figure_index <- if (anyNA(figure_indices)) {
      NA_integer_
    } else {
      max(figure_indices)
    }
    if (
      is.na(marker_index) ||
        is.na(de_index) ||
        is.na(last_figure_index) ||
        marker_index <= last_figure_index ||
        de_index <= marker_index
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": marker/DE stage order mismatch; expected mg-markers then mg-de after all figure stages; got ",
          paste(observed_stages, collapse = ", ")
        )
      )
    }
    terminal_stage_contracts <- list(
      "render-notebook" = list(
        command = "command: 'quarto' 'render' 'notebook/sc_analysis.qmd'",
        expects = "notebook/sc_analysis.html"
      ),
      "tripwires" = list(
        command = "command: 'Rscript' 'tools/run-tripwires.R'",
        expects = "tools/run-tripwires.R"
      )
    )
    if (
      !identical(
        tail(observed_stages, 3L),
        c("mg-de", "render-notebook", "tripwires")
      )
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": terminal stage order mismatch; expected mg-de then render-notebook then tripwires; got ",
          paste(tail(observed_stages, 3L), collapse = ", ")
        )
      )
    }
    for (terminal_stage in names(terminal_stage_contracts)) {
      contract <- terminal_stage_contracts[[terminal_stage]]
      record <- stage_records[[terminal_stage]]
      if (!identical(record$command, contract$command)) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            terminal_stage,
            paste0("exact command ", shQuote(contract$command)),
            record$command
          )
        )
      }
      output_paths <- stage_output_paths(record)
      if (
        length(output_paths) != 1L ||
          !path_has_suffix(output_paths[[1L]], contract$expects)
      ) {
        problems <- c(
          problems,
          stage_problem(
            invocation$label,
            terminal_stage,
            paste0("exact expected output path ", shQuote(contract$expects)),
            record$expects
          )
        )
      }
    }
  }

  for (invocation in failed_invocations) {
    result <- run_pipeline(invocation$args)
    if (identical(result$status, 0L)) {
      problems <- c(
        problems,
        paste0(invocation$label, ": expected a nonzero exit status, got 0")
      )
    }
    has_diagnostic <- invocation$expected %in%
      result$output ||
      paste0("Error: ", invocation$expected) %in% result$output
    if (!has_diagnostic) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": diagnostic mismatch; expected ",
          shQuote(invocation$expected),
          "; got ",
          format_output(result$output)
        )
      )
    }
  }

  missing_input <- tempfile(
    pattern = "espi-tripwire-missing-explicit-input-",
    tmpdir = tempdir(),
    fileext = ".rds"
  )
  missing_input_expected <- paste0(
    "Pipeline input(s) do not exist: ",
    missing_input
  )
  missing_input_result <- run_pipeline(c("--input", missing_input))
  if (identical(missing_input_result$status, 0L)) {
    problems <- c(
      problems,
      "missing explicit input no-side-effect: expected a nonzero exit status, got 0"
    )
  }
  if (
    !identical(
      missing_input_result$output,
      c(paste0("Error: ", missing_input_expected), "Execution halted")
    )
  ) {
    problems <- c(
      problems,
      paste0(
        "missing explicit input no-side-effect: diagnostic mismatch; expected ",
        format_output(c(
          paste0("Error: ", missing_input_expected),
          "Execution halted"
        )),
        "; got ",
        format_output(missing_input_result$output)
      )
    )
  }
  if (file.exists(missing_input)) {
    problems <- c(
      problems,
      paste0(
        "missing explicit input no-side-effect: expected path to remain absent: ",
        shQuote(missing_input)
      )
    )
  }

  if (length(problems) > 0L) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }
  pass(
    slug,
    "Public pipeline dry-run exposes the complete deterministic plan and execution preflight has no side effects."
  )
}

# Operational boundary: the public just recipes must be a thin, observable
# interface over the same deterministic pipeline runner as the direct CLI.
# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_just_public_interface <- function(root) {
  slug <- "just-public-interface"
  just <- Sys.which("just")
  rscript <- Sys.which("Rscript")
  if (identical(unname(just), "")) {
    return(fail(
      slug,
      "just is unavailable, so the public recipe interface cannot be executed."
    ))
  }
  if (identical(unname(rscript), "")) {
    return(fail(
      slug,
      "Rscript is unavailable, so the public recipe interface cannot be executed."
    ))
  }

  script <- file.path(root, "scripts", "run-pipeline.R")
  if (!file.exists(script)) {
    return(fail(slug, "scripts/run-pipeline.R is missing."))
  }

  data_root_env <- Sys.getenv("MEGAN_SC_DATA_DIR", unset = "")
  box_path_env <- Sys.getenv("BOX_PATH", unset = "")
  data_roots <- unique(c(
    data_root_env,
    if (nzchar(box_path_env)) file.path(box_path_env, "megan_sc_data") else "",
    box_path_env,
    file.path(path.expand("~/Library/CloudStorage/Box-Box"), "megan_sc_data"),
    file.path(root, "data")
  ))
  # ANALYSIS_OK[filtering]: discard empty DATA_ROOT candidates before constructing input directories; the resulting candidate set is checked by the missing-input tripwire.
  data_roots <- data_roots[nzchar(data_roots)]
  input_dirs <- unique(c(
    file.path(data_roots, "seurat_objects", "input"),
    file.path(data_roots, "input")
  ))
  explicit_candidates <- unique(unlist(lapply(
    input_dirs,
    function(path) {
      if (!dir.exists(path)) {
        return(character())
      }
      list.files(path, pattern = "\\.rds$", full.names = TRUE)
    }
  )))
  # ANALYSIS_OK[filtering]: retain only existing explicit RDS candidates so the missing-input contract can assert a controlled failure when none remain.
  explicit_candidates <- explicit_candidates[file.exists(explicit_candidates)]
  if (length(explicit_candidates) == 0L) {
    return(fail(
      slug,
      "No existing explicit RDS input was found for the public recipe exercise."
    ))
  }
  explicit_input <- sort(explicit_candidates)[[1L]]

  run_process <- function(program, args) {
    output <- tryCatch(
      # ANALYSIS_OK[warning-suppression]: command probes intentionally capture expected subprocess warnings as output/status; each caller asserts the failure contract.
      suppressWarnings(system2(
        program,
        args,
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(error) structure(conditionMessage(error), status = 1L)
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }
    list(status = as.integer(status), output = unname(as.character(output)))
  }

  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  plan_lines <- function(output) {
    output[grepl(
      paste0(
        "^(mode: dry-run|input_source: |overwrite: |",
        "source_cluster_column: |mg_cluster_column: |mg_pca_dims: |",
        "first_stage: |final_stage: |stage_count: |stage: |",
        "command: |expects: )"
      ),
      output,
      perl = TRUE
    )]
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  stage_names <- function(plan) {
    sub("^stage: ", "", plan[startsWith(plan, "stage: ")])
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  command_label <- function(label, side) {
    paste(label, side, sep = " / ")
  }

  successful_invocations <- list(
    list(
      label = "counts-qc default",
      just_args = c("--quiet", "run-dry-run"),
      direct_args = c(shQuote(script), "--dry-run"),
      input_source = "counts-qc",
      overwrite = FALSE,
      expected_stage_count = 24L
    ),
    list(
      label = "legacy",
      just_args = c("--quiet", "run-dry-run", "legacy", "false"),
      direct_args = c(shQuote(script), "--dry-run", "--input-source", "legacy"),
      input_source = "legacy",
      overwrite = FALSE,
      expected_stage_count = 22L
    ),
    list(
      label = "explicit existing RDS",
      just_args = c(
        "--quiet",
        "run-dry-run",
        shQuote(explicit_input),
        "false"
      ),
      direct_args = c(
        shQuote(script),
        "--dry-run",
        "--input",
        shQuote(explicit_input)
      ),
      input_source = "explicit",
      overwrite = FALSE,
      expected_stage_count = 22L
    ),
    list(
      label = "counts-qc overwrite",
      just_args = c("--quiet", "run-dry-run", "counts-qc", "true"),
      direct_args = c(shQuote(script), "--dry-run", "--overwrite"),
      input_source = "counts-qc",
      overwrite = TRUE,
      expected_stage_count = 24L
    )
  )

  problems <- character()
  for (invocation in successful_invocations) {
    direct <- run_process(rscript, invocation$direct_args)
    recipe <- run_process(just, invocation$just_args)
    if (!identical(direct$status, 0L)) {
      problems <- c(
        problems,
        sprintf(
          "%s: direct CLI exited %d: %s",
          command_label(invocation$label, "direct"),
          direct$status,
          paste(direct$output, collapse = " | ")
        )
      )
      next
    }
    if (!identical(recipe$status, 0L)) {
      problems <- c(
        problems,
        sprintf(
          "%s: recipe exited %d: %s",
          command_label(invocation$label, "just"),
          recipe$status,
          paste(recipe$output, collapse = " | ")
        )
      )
      next
    }

    direct_plan <- plan_lines(direct$output)
    recipe_plan <- plan_lines(recipe$output)
    if (!identical(recipe_plan, direct_plan)) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": just plan differs from direct CLI contract; direct=",
          paste(shQuote(direct_plan), collapse = ", "),
          "; just=",
          paste(shQuote(recipe_plan), collapse = ", ")
        )
      )
      next
    }

    expected_header <- c(
      "mode: dry-run",
      paste0("input_source: ", invocation$input_source),
      paste0("overwrite: ", tolower(as.character(invocation$overwrite)))
    )
    if (
      length(direct_plan) < length(expected_header) ||
        !identical(direct_plan[seq_along(expected_header)], expected_header)
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": direct CLI header does not match the expected source/overwrite contract."
        )
      )
      next
    }

    observed_stages <- stage_names(direct_plan)
    observed_count <- length(observed_stages)
    if (!identical(observed_count, invocation$expected_stage_count)) {
      problems <- c(
        problems,
        sprintf(
          "%s: expected %d stages, got %d",
          invocation$label,
          invocation$expected_stage_count,
          observed_count
        )
      )
    }
    if (
      length(observed_stages) == 0L ||
        !identical(
          observed_stages[[1L]],
          if (identical(invocation$input_source, "counts-qc")) {
            "process-counts"
          } else {
            "preprocess-source"
          }
        )
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": source-dependent first stage is wrong: ",
          paste(observed_stages, collapse = ", ")
        )
      )
    }
    has_counts_stages <- all(
      c("process-counts", "qc-filtering") %in% observed_stages
    )
    if (
      identical(invocation$input_source, "counts-qc") &&
        !isTRUE(has_counts_stages)
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": counts-qc plan is missing process-counts/qc-filtering stages."
        )
      )
    }
    if (
      !identical(invocation$input_source, "counts-qc") &&
        isTRUE(has_counts_stages)
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": legacy/explicit plan unexpectedly includes counts-qc stages."
        )
      )
    }
    if (
      length(observed_stages) < MIN_EXPECTED_FINAL_STAGES ||
        !identical(
          tail(observed_stages, 3L),
          c("mg-de", "render-notebook", "tripwires")
        )
    ) {
      problems <- c(
        problems,
        paste0(
          invocation$label,
          ": final stages are not mg-de/render-notebook/tripwires: ",
          paste(tail(observed_stages, 3L), collapse = ", ")
        )
      )
    }
  }

  invalid <- run_process(
    just,
    c("--quiet", "run-dry-run", "counts-qc", "not-a-boolean")
  )
  if (identical(invalid$status, 0L)) {
    problems <- c(
      problems,
      "invalid overwrite value was accepted by just run-dry-run."
    )
  } else if (
    !any(grepl("overwrite", invalid$output, ignore.case = TRUE)) ||
      !any(grepl("true|false", invalid$output, ignore.case = TRUE))
  ) {
    problems <- c(
      problems,
      paste0(
        "invalid overwrite value failed without an overwrite boolean diagnostic: ",
        paste(invalid$output, collapse = " | ")
      )
    )
  }

  listing <- run_process(just, c("--list"))
  if (!identical(listing$status, 0L)) {
    problems <- c(
      problems,
      paste0(
        "just --list exited ",
        listing$status,
        ": ",
        paste(listing$output, collapse = " | ")
      )
    )
  } else {
    recipe_lines <- listing$output[grepl(
      "^[[:space:]]+[A-Za-z0-9_-]+([[:space:]]|$)",
      listing$output,
      perl = TRUE
    )]
    recipe_names <- sub(
      "^[[:space:]]+([A-Za-z0-9_-]+).*",
      "\\1",
      recipe_lines
    )
    required_names <- c("run", "run-dry-run", "preprocess", "preprocess-one")
    missing_names <- setdiff(required_names, recipe_names)
    if (length(missing_names) > 0L) {
      problems <- c(
        problems,
        paste0(
          "just --list is missing recipe(s): ",
          paste(missing_names, collapse = ", ")
        )
      )
    } else {
      low_level_index <- min(match(
        c("preprocess", "preprocess-one"),
        recipe_names
      ))
      if (
        match("run", recipe_names) >= low_level_index ||
          match("run-dry-run", recipe_names) >= low_level_index
      ) {
        problems <- c(
          problems,
          "just --list places run/run-dry-run after low-level preprocessing recipes."
        )
      }
    }
  }

  if (length(problems) > 0L) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }
  pass(
    slug,
    "Public just run/run-dry-run recipes match direct CLI dry-run plans for counts-qc, legacy, and an existing explicit RDS, including overwrite validation and recipe ordering."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_report_values_freshness <- function(root) {
  # Scientific boundary: the rendered report must not be older than the source
  # prose or figures that define the claimed HVG and DimHeatmap parameters.
  slug <- "report-values-freshness"
  qmd <- file.path(root, "notebook", "sc_analysis.qmd")
  html <- file.path(root, "notebook", "sc_analysis.html")
  if (!file.exists(qmd) || !file.exists(html)) {
    missing <- c(qmd, html)[!file.exists(c(qmd, html))]
    return(fail(
      slug,
      paste(
        "Missing report artifact(s):",
        paste(basename(missing), collapse = ", ")
      )
    ))
  }

  qmd_lines <- read_text(qmd)
  qmd_text <- squash(qmd_lines)

  figure_refs <- unique(unlist(regmatches(
    qmd_text,
    gregexpr("figures/[^)\"'{}[:space:]]+\\.png", qmd_text, perl = TRUE)
  )))
  figure_paths <- file.path(dirname(qmd), figure_refs)
  missing_figure_refs <- figure_refs[!file.exists(figure_paths)]
  missing_figures <- file.path(dirname(qmd), missing_figure_refs)

  problems <- character()
  if (length(missing_figures) > 0) {
    problems <- c(
      problems,
      paste(
        "missing QMD figure reference(s):",
        paste(file.path("notebook", missing_figure_refs), collapse = ", ")
      )
    )
  }

  html_mtime <- mtime(html)
  qmd_mtime <- mtime(qmd)
  if (is.na(html_mtime) || is.na(qmd_mtime) || html_mtime <= qmd_mtime) {
    problems <- c(
      problems,
      "notebook/sc_analysis.html is not newer than notebook/sc_analysis.qmd"
    )
  }

  existing_figures <- figure_paths[file.exists(figure_paths)]
  target_paths <- vapply(existing_figures, resolve_figure_target, character(1))
  missing_targets <- target_paths[!file.exists(target_paths)]
  if (length(missing_targets) > 0) {
    problems <- c(
      problems,
      paste(
        "missing figure symlink target(s):",
        paste(missing_targets, collapse = ", ")
      )
    )
  }

  existing_targets <- target_paths[file.exists(target_paths)]
  stale_targets <- existing_targets[
    !is.na(mtime(existing_targets)) & html_mtime <= mtime(existing_targets)
  ]
  if (length(stale_targets) > 0) {
    rel <- sub(paste0("^", root, "/"), "", stale_targets)
    problems <- c(
      problems,
      paste(
        "notebook/sc_analysis.html is not newer than figure target(s):",
        paste(rel, collapse = ", ")
      )
    )
  }

  mentions_top20 <- grepl(
    "top[[:space:]]+20[[:space:]]+(retained[[:space:]]+)?HVGs?",
    qmd_text,
    ignore.case = TRUE,
    perl = TRUE
  )
  mentions_500_cells <- grepl(
    "500[[:space:]]+cells",
    qmd_text,
    ignore.case = TRUE,
    perl = TRUE
  )
  if (!mentions_top20) {
    problems <- c(
      problems,
      "QMD prose does not mention the top 20 HVGs labeling contract"
    )
  }
  if (!mentions_500_cells) {
    problems <- c(
      problems,
      "QMD prose does not mention the 500 cells DimHeatmap contract"
    )
  }

  if (length(problems) > 0) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }
  pass(
    slug,
    sprintf(
      "HTML is newer than QMD and %d referenced figure target(s); prose states top 20 HVGs and 500 cells.",
      length(existing_targets)
    )
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_missing_counts_file <- function(root) {
  # Scientific boundary: an explicit missing input must abort instead of silently
  # falling back to a default object or stale cache.
  slug <- "missing-counts-file"
  script <- file.path(root, "scripts", "03-preprocess.R")
  if (!file.exists(script)) {
    return(fail(slug, "scripts/03-preprocess.R is missing."))
  }
  rscript <- Sys.which("Rscript")
  if (identical(unname(rscript), "")) {
    return(fail(
      slug,
      "Rscript is unavailable, so the missing-input behavior cannot be executed."
    ))
  }

  missing_input <- file.path(
    tempdir(),
    sprintf("espi-missing-input-%s.rds", Sys.getpid())
  )
  if (file.exists(missing_input)) {
    unlink(missing_input)
  }

  output <- tryCatch(
    system2(
      rscript,
      c(script, "--input", missing_input, "--normalization", "log1p"),
      stdout = TRUE,
      stderr = TRUE
    ),
    warning = function(w) structure(conditionMessage(w), status = 1L),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  text <- paste(output, collapse = "\n")

  if (identical(as.integer(status), 0L)) {
    return(fail(
      slug,
      "preprocess-sobj.R exited 0 for a deliberately missing --input path."
    ))
  }
  if (
    !grepl(missing_input, text, fixed = TRUE) &&
      !grepl(
        "No such file|cannot open|readRDS",
        text,
        ignore.case = TRUE,
        perl = TRUE
      )
  ) {
    msg <- paste(
      "preprocess-sobj.R failed non-zero, but not at the missing input boundary. First output:",
      substr(gsub("[\r\n]+", " ", text), 1, 240)
    )
    return(fail(slug, msg))
  }

  pass(
    slug,
    "preprocess-sobj.R returns non-zero for a deliberately missing --input path and surfaces the missing-file boundary."
  )
}
# Execute script 01 against read-only count-directory links and a metadata
# fixture with one sample removed. The explicit output path keeps Box artifacts
# out of this fault injection.
# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_missing_metadata_sample <- function(root) {
  slug <- "missing-metadata-sample"
  script <- file.path(root, "scripts", "01-process-counts.R")
  rscript <- Sys.which("Rscript")
  if (!file.exists(script)) {
    return(fail(slug, "scripts/01-process-counts.R is missing."))
  }
  if (identical(unname(rscript), "")) {
    return(fail(
      slug,
      "Rscript is unavailable, so metadata reconciliation cannot be executed."
    ))
  }

  data_root <- Sys.getenv("DATA_ROOT_DIR", unset = "")
  if (!nzchar(data_root)) {
    data_root <- file.path(
      path.expand("~/Library/CloudStorage/Box-Box"),
      "megan_sc_data"
    )
  }
  raw_dir <- file.path(data_root, "data", "input", "Raw Matrices")
  metadata_path <- file.path(raw_dir, "Sample_Metadata_MS1.txt")
  if (!dir.exists(raw_dir) || !file.exists(metadata_path)) {
    return(skip(
      slug,
      "Scratch reconciliation requires the production six-sample Raw Matrices directory and Sample_Metadata_MS1.txt."
    ))
  }

  metadata <- tryCatch(
    utils::read.delim(
      metadata_path,
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    # ANALYSIS_OK[optional-input]: missing metadata is an intentional fault-injection input; this branch records the controlled probe outcome below.
    error = function(e) {
      message(sprintf("Metadata fixture read failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(metadata) || nrow(metadata) < MIN_METADATA_ROWS) {
    return(skip(
      slug,
      "Production metadata is unavailable or does not contain at least two Sample rows."
    ))
  }
  samples <- as.character(metadata$Sample)
  sample_dirs <- file.path(raw_dir, samples)
  if (any(!dir.exists(sample_dirs))) {
    return(skip(
      slug,
      "Production count directories do not reconcile with metadata, so the fault injection has no valid baseline."
    ))
  }

  scratch <- tempfile("espi-missing-metadata-")
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)
  scratch_raw <- file.path(scratch, "Raw Matrices")
  dir.create(scratch_raw, recursive = TRUE, showWarnings = FALSE)
  for (sample in samples) {
    file.symlink(
      normalizePath(
        file.path(raw_dir, sample),
        winslash = "/",
        mustWork = TRUE
      ),
      file.path(scratch_raw, sample)
    )
  }
  missing_sample <- samples[[1L]]
  utils::write.table(
    metadata[-1L, , drop = FALSE],
    file.path(scratch_raw, "Sample_Metadata_MS1.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  checkpoint_log <- file.path(scratch, "checkpoints.tsv")
  drop_ledger <- file.path(scratch, "drops.tsv")
  output_path <- file.path(scratch, "sobj_raw.rds")
  env_names <- c("CHECKPOINT_LOG", "DROP_LEDGER", "STOP_AFTER_CHECKPOINT")
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  on.exit(
    for (i in seq_along(env_names)) {
      if (is.na(old_env[[i]])) {
        Sys.unsetenv(env_names[[i]])
      } else {
        Sys.setenv(structure(old_env[[i]], names = env_names[[i]]))
      }
    },
    add = TRUE
  )
  Sys.setenv(
    CHECKPOINT_LOG = checkpoint_log,
    DROP_LEDGER = drop_ledger,
    STOP_AFTER_CHECKPOINT = "samples_reconciled"
  )
  output <- tryCatch(
    system2(
      rscript,
      c(
        script,
        "--raw-counts-dir",
        shQuote(scratch_raw),
        "--metadata",
        shQuote(file.path(scratch_raw, "Sample_Metadata_MS1.txt")),
        "--output",
        shQuote(output_path)
      ),
      stdout = TRUE,
      stderr = TRUE,
      env = c(
        paste0("CHECKPOINT_LOG=", checkpoint_log),
        paste0("DROP_LEDGER=", drop_ledger),
        "STOP_AFTER_CHECKPOINT=samples_reconciled"
      )
    ),
    warning = function(w) structure(conditionMessage(w), status = 1L),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  checkpoint_text <- if (file.exists(checkpoint_log)) {
    read_text(checkpoint_log)
  } else {
    character()
  }
  has_reconciled <- any(grepl(
    "\tsamples_reconciled\t",
    checkpoint_text,
    fixed = TRUE
  ))
  ledger_ok <- FALSE
  if (file.exists(drop_ledger)) {
    ledger <- tryCatch(
      utils::read.delim(
        drop_ledger,
        sep = "\t",
        quote = "\"",
        stringsAsFactors = FALSE
      ),
      # ANALYSIS_OK[optional-input]: malformed reconciliation ledger is an intentional fault-injection input; the tripwire verifies the resulting failure record.
      error = function(e) {
        message(sprintf(
          "Reconciliation ledger read failed: %s",
          conditionMessage(e)
        ))
        NULL
      }
    )
    if (
      !is.null(ledger) &&
        all(c("sample_id", "reason", "allowed_by_policy") %in% names(ledger))
    ) {
      allowed <- tolower(trimws(as.character(ledger$allowed_by_policy)))
      ledger_ok <- any(allowed %in% c("false", "0", "no"))
    }
  }
  if (!ledger_ok && file.exists(drop_ledger)) {
    ledger_ok <- length(read_text(drop_ledger)) > 1L &&
      any(grepl(
        "FALSE",
        read_text(drop_ledger),
        ignore.case = TRUE,
        fixed = TRUE
      ))
  }
  if (!identical(as.integer(status), 0L) && !has_reconciled && ledger_ok) {
    return(pass(
      slug,
      sprintf(
        "Missing metadata sample %s fails reconciliation, emits no samples_reconciled checkpoint, and records disallowed ledger evidence.",
        missing_sample
      )
    ))
  }
  fail(
    slug,
    sprintf(
      "Expected non-zero reconciliation failure with no samples_reconciled checkpoint and disallowed ledger evidence (exit=%s, checkpoint=%s, ledger=%s, ledger_file=%s).",
      status,
      has_reconciled,
      ledger_ok,
      file.exists(drop_ledger)
    )
  )
}


# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_heatmap_missing_input <- function(root) {
  # Operational boundary: heatmap plotting must fail at the explicit missing
  # input instead of falling back to a default Seurat object or writing partial
  # figure/table/notebook artifacts.
  slug <- "heatmap-missing-input"
  script <- file.path(root, "scripts", "10-plot-cluster-marker-heatmaps.R")
  if (!file.exists(script)) {
    return(fail(slug, "scripts/10-plot-cluster-marker-heatmaps.R is missing."))
  }
  rscript <- Sys.which("Rscript")
  if (identical(unname(rscript), "")) {
    return(fail(
      slug,
      "Rscript is unavailable, so the heatmap missing-input behavior cannot be executed."
    ))
  }

  # ANALYSIS_OK[file-freshness-tripwire]: mtime snapshots are the explicit side-effect signal; the before/after comparison below detects unexpected production writes.
  snapshot_dir_state <- function(dir) {
    if (!dir.exists(dir)) {
      return(character())
    }
    files <- list.files(
      dir,
      all.files = TRUE,
      no.. = TRUE,
      recursive = TRUE,
      full.names = TRUE
    )
    files <- sort(normalizePath(files, winslash = "/", mustWork = FALSE))
    if (length(files) == 0L) {
      return(character())
    }
    info <- file.info(files)
    state <- paste(info$size, as.numeric(info$mtime), sep = ":")
    names(state) <- files
    state
  }

  scratch_root <- tempfile("espi-heatmap-missing-")
  scratch_out_dir <- file.path(scratch_root, "figures")
  dir.create(scratch_root, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(scratch_root, recursive = TRUE, force = TRUE), add = TRUE)
  missing_input <- file.path(scratch_root, "missing-cluster-object.rds")
  table_dir <- analysis_table_annotation_dir(root)
  notebook_figure_dir <- file.path(root, "notebook", "figures")
  before_tables <- snapshot_dir_state(table_dir)
  before_notebook <- snapshot_dir_state(notebook_figure_dir)

  output <- tryCatch(
    # ANALYSIS_OK[warning-suppression]: missing-input subprocess probes intentionally capture warnings as status/output; the tripwire asserts both failure and scratch cleanup.
    suppressWarnings(system2(
      rscript,
      c(
        script,
        "--input",
        missing_input,
        "--dims",
        "999",
        "--resolution",
        "tripwire_missing_input",
        "--n-perm",
        "1",
        "--out-dir",
        scratch_out_dir
      ),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  text <- paste(output, collapse = "\n")

  after_tables <- snapshot_dir_state(table_dir)
  after_notebook <- snapshot_dir_state(notebook_figure_dir)
  scratch_files <- if (dir.exists(scratch_out_dir)) {
    list.files(
      scratch_out_dir,
      all.files = TRUE,
      no.. = TRUE,
      recursive = TRUE,
      full.names = TRUE
    )
  } else {
    character()
  }
  table_changes <- union(
    setdiff(names(after_tables), names(before_tables)),
    names(after_tables)[
      names(after_tables) %in%
        names(before_tables) &
        after_tables != before_tables[names(after_tables)]
    ]
  )
  notebook_changes <- union(
    setdiff(names(after_notebook), names(before_notebook)),
    names(after_notebook)[
      names(after_notebook) %in%
        names(before_notebook) &
        after_notebook != before_notebook[names(after_notebook)]
    ]
  )
  notebook_links <- list.files(
    notebook_figure_dir,
    all.files = TRUE,
    no.. = TRUE,
    full.names = TRUE
  )
  link_targets <- Sys.readlink(notebook_links)
  resolved_link_targets <- ifelse(
    grepl("^/", link_targets),
    link_targets,
    file.path(dirname(notebook_links), link_targets)
  )
  scratch_links <- notebook_links[
    nzchar(link_targets) &
      grepl(
        scratch_root,
        normalizePath(
          resolved_link_targets,
          winslash = "/",
          mustWork = FALSE
        ),
        fixed = TRUE
      )
  ]

  problems <- character()
  if (identical(as.integer(status), 0L)) {
    problems <- c(
      problems,
      "plot-cluster-marker-heatmaps.R exited 0 for a deliberately missing --input path"
    )
  }
  if (!grepl("Input Seurat object does not exist", text, fixed = TRUE)) {
    problems <- c(
      problems,
      paste(
        "missing expected heatmap input error. First output:",
        substr(gsub("[\r\n]+", " ", text), 1, 240)
      )
    )
  }
  if (length(scratch_files) > 0L) {
    problems <- c(
      problems,
      paste(
        "scratch heatmap output(s) were created:",
        paste(basename(scratch_files), collapse = ", ")
      )
    )
  }
  if (length(table_changes) > 0L) {
    problems <- c(
      problems,
      paste(
        "TABLE_DIR/annotation changed during missing-input run:",
        paste(basename(table_changes), collapse = ", ")
      )
    )
  }
  if (length(notebook_changes) > 0L) {
    problems <- c(
      problems,
      paste(
        "notebook figure link(s) changed during missing-input run:",
        paste(basename(notebook_changes), collapse = ", ")
      )
    )
  }
  if (length(scratch_links) > 0L) {
    problems <- c(
      problems,
      paste(
        "notebook figure symlink(s) point into scratch output:",
        paste(basename(scratch_links), collapse = ", ")
      )
    )
  }

  if (length(problems) > 0L) {
    return(fail(slug, compact_problem_list(problems)))
  }

  pass(
    slug,
    "plot-cluster-marker-heatmaps.R fails non-zero at the missing --input boundary and leaves scratch figures, fixed annotation tables, and notebook figure links unchanged."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_p27_rng_state_preservation <- function(root) {
  # Statistical helper boundary: p27 enrichment permutations may be seeded for
  # determinism, but must not advance or replace the caller's RNG state.
  slug <- "rng-state-preservation"
  required_packages <- c("devtools", "Seurat", "SeuratObject")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    return(fail(
      slug,
      paste(
        "Missing required package(s) for RNG tripwire:",
        paste(missing_packages, collapse = ", ")
      )
    ))
  }

  random_seed_existed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  original_random_seed <- if (random_seed_existed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(
    restore_random_seed(original_random_seed, random_seed_existed),
    add = TRUE
  )

  devtools::load_all(root, export_all = FALSE, quiet = TRUE)
  compute_cluster_p27_enrichment <- get(
    "compute_cluster_p27_enrichment",
    envir = asNamespace("ESPI"),
    inherits = FALSE
  )
  sobj <- make_tripwire_p27_sobj()

  # ANALYSIS_OK[random-seed-only]: deterministic fixture seed is required to compare p27 enrichment outputs; restore_random_seed() verifies caller RNG preservation.
  set.seed(20260705)
  expected_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  result_one <- compute_cluster_p27_enrichment(
    sobj,
    "cluster",
    layer = "pflog",
    condition_col = "Condition",
    n_perm = 80L,
    seed = 4242L
  )
  actual_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (!identical(actual_seed, expected_seed)) {
    return(fail(
      slug,
      "compute_cluster_p27_enrichment() changed .Random.seed for an existing caller seed."
    ))
  }

  # ANALYSIS_OK[random-seed-only]: second deterministic fixture seed is required to verify repeatability; restore_random_seed() verifies caller RNG preservation.
  set.seed(20260706)
  result_two <- compute_cluster_p27_enrichment(
    sobj,
    "cluster",
    layer = "pflog",
    condition_col = "Condition",
    n_perm = 80L,
    seed = 4242L
  )
  if (
    !isTRUE(all.equal(
      result_one$z_score,
      result_two$z_score,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))
  ) {
    return(fail(
      slug,
      "compute_cluster_p27_enrichment() returned non-deterministic z-scores for the same seed."
    ))
  }
  if (anyNA(result_one$z_score) || any(!is.finite(result_one$z_score))) {
    return(fail(slug, "Synthetic p27 enrichment returned non-finite z-scores."))
  }
  enriched <- result_one$z_score[match("1", result_one$cluster)]
  depleted <- result_one$z_score[match("2", result_one$cluster)]
  if (!isTRUE(enriched > 0 && depleted < 0)) {
    return(fail(
      slug,
      "Synthetic elevated/depleted p27 clusters did not return positive/negative z-scores."
    ))
  }

  pass(
    slug,
    "compute_cluster_p27_enrichment() preserves caller RNG state exactly and returns deterministic positive/negative synthetic p27 z-scores."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_mycelium_provenance_semantics <- function(root) {
  slug <- "mycelium-provenance-semantics"
  registry <- file.path(root, ".living", "log", "LOG_REGISTRY.md")

  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  validate_registry <- function(
    registry_path,
    session_ids,
    strict_ids = character()
  ) {
    if (!file.exists(registry_path)) {
      return(sprintf("%s is missing", basename(registry_path)))
    }
    registry_lines <- read_text(registry_path)
    table_lines <- grep("^\\|", registry_lines)
    if (length(table_lines) < MIN_MARKDOWN_TABLE_LINES) {
      return(sprintf(
        "%s does not contain a markdown table",
        basename(registry_path)
      ))
    }
    header <- parse_markdown_table_row(registry_lines[[table_lines[[1L]]]])
    required <- c(
      "Session ID",
      "Summary",
      "Key Outputs",
      "Status",
      "Tags",
      "Log",
      "Duration",
      "Files Changed"
    )
    missing_columns <- setdiff(required, header)
    if (length(missing_columns) > 0L) {
      return(paste(
        "registry is missing required column(s):",
        paste(missing_columns, collapse = ", ")
      ))
    }
    rows <- list()
    problems <- character()
    for (line_no in table_lines[-REGISTRY_HEADER_ROW_INDEX]) {
      fields <- parse_markdown_table_row(registry_lines[[line_no]])
      if (length(fields) != length(header)) {
        next
      }
      row <- stats::setNames(fields, header)
      if (row[["Session ID"]] %in% session_ids) {
        rows[[row[["Session ID"]]]] <- list(row = row, line_no = line_no)
      }
    }
    missing_ids <- setdiff(session_ids, names(rows))
    if (length(missing_ids) > 0L) {
      problems <- c(
        problems,
        paste("missing registry row(s):", paste(missing_ids, collapse = ", "))
      )
    }
    repo_root <- dirname(dirname(dirname(registry_path)))
    for (session_id in intersect(session_ids, names(rows))) {
      entry <- rows[[session_id]]
      row <- entry$row
      line_no <- entry$line_no
      strict <- session_id %in% strict_ids
      if (!identical(tolower(trimws(row[["Status"]])), "complete")) {
        problems <- c(
          problems,
          sprintf("L%d %s status is not complete", line_no, session_id)
        )
      }
      for (field in c("Summary", "Key Outputs", "Tags")) {
        value <- trimws(row[[field]])
        if (!nzchar(value) || identical(value, "—")) {
          problems <- c(
            problems,
            sprintf("L%d %s has empty %s", line_no, session_id, field)
          )
        }
      }
      if (is_file_list_only_summary(row[["Summary"]])) {
        problems <- c(
          problems,
          sprintf("L%d %s summary is file-list-only", line_no, session_id)
        )
      }
      log_target <- extract_markdown_link_target(row[["Log"]])
      if (is.na(log_target) || !nzchar(log_target)) {
        problems <- c(
          problems,
          sprintf("L%d %s lacks a linked log", line_no, session_id)
        )
        next
      }
      log_path <- file.path(dirname(registry_path), log_target)
      if (!file.exists(log_path)) {
        problems <- c(
          problems,
          sprintf(
            "L%d %s linked log is missing: %s",
            line_no,
            session_id,
            log_target
          )
        )
        next
      }
      log_lines <- read_text(log_path)
      for (field in c("ended", "duration_minutes", "files_changed")) {
        value <- frontmatter_scalar(log_lines, field)
        if (is.na(value) || !nzchar(value) || identical(value, "NA")) {
          problems <- c(
            problems,
            sprintf(
              "L%d %s linked log lacks populated frontmatter %s",
              line_no,
              session_id,
              field
            )
          )
        }
      }
      duration_match <- regexec(
        "([0-9]+)[[:space:]]*m",
        row[["Duration"]],
        perl = TRUE
      )
      duration_hit <- regmatches(row[["Duration"]], duration_match)[[1L]]
      # ANALYSIS_OK[warning-suppression]: malformed optional frontmatter is parsed as NA for a validation comparison; the surrounding check reports the discrepancy.
      front_duration <- suppressWarnings(as.numeric(frontmatter_scalar(
        log_lines,
        "duration_minutes"
      )))
      if (
        length(duration_hit) == REGEX_CAPTURE_COUNT &&
          is.finite(front_duration) &&
          !identical(
            as.numeric(duration_hit[[REGEX_CAPTURE_VALUE]]),
            front_duration
          )
      ) {
        problems <- c(
          problems,
          sprintf(
            "L%d %s duration disagrees with linked log",
            line_no,
            session_id
          )
        )
      }
      files_match <- regexec(
        "^[[:space:]]*([0-9]+)",
        row[["Files Changed"]],
        perl = TRUE
      )
      files_hit <- regmatches(row[["Files Changed"]], files_match)[[1L]]
      # ANALYSIS_OK[warning-suppression]: malformed optional frontmatter is parsed as NA for a validation comparison; the surrounding check reports the discrepancy.
      front_files <- suppressWarnings(as.numeric(frontmatter_scalar(
        log_lines,
        "files_changed"
      )))
      if (
        length(files_hit) == REGEX_CAPTURE_COUNT &&
          is.finite(front_files) &&
          !identical(as.numeric(files_hit[[REGEX_CAPTURE_VALUE]]), front_files)
      ) {
        problems <- c(
          problems,
          sprintf(
            "L%d %s files_changed disagrees with linked log",
            line_no,
            session_id
          )
        )
      }
      if (!strict) {
        next
      }
      for (section in c("Session Summary", "Key Outputs", "Status")) {
        if (
          !any(grepl(
            sprintf("^##[[:space:]]+%s[[:space:]]*$", section),
            log_lines,
            perl = TRUE
          ))
        ) {
          problems <- c(
            problems,
            sprintf(
              "L%d %s linked log lacks ## %s",
              line_no,
              session_id,
              section
            )
          )
        }
      }
      files_start <- grep(
        "^##[[:space:]]+Files Modified[[:space:]]*$",
        log_lines,
        perl = TRUE
      )
      if (length(files_start) > 0L) {
        next_heading <- which(
          seq_along(log_lines) > files_start[[1L]] &
            grepl("^#{1,6}[[:space:]]+", log_lines, perl = TRUE)
        )
        files_end <- if (length(next_heading) > 0L) {
          next_heading[[1L]] - 1L
        } else {
          length(log_lines)
        }
        file_lines <- log_lines[seq.int(files_start[[1L]] + 1L, files_end)]
        for (file_line in file_lines[grepl("^[-*][[:space:]]+", file_lines)]) {
          target <- trimws(sub("^[-*][[:space:]]+", "", file_line))
          target <- sub("^`(.*)`$", "\\1", target)
          deleted <- grepl(
            "\\b(deleted|removed)\\b",
            target,
            ignore.case = TRUE,
            perl = TRUE
          )
          target <- sub("^[Dd]eleted:[[:space:]]*", "", target)
          if (grepl("^(local://|\\.omp/|/Users/.*/\\.omp/)", target)) {
            next
          }
          target_path <- if (grepl("^/", target)) {
            target
          } else {
            file.path(repo_root, target)
          }
          if (!deleted && !file.exists(target_path)) {
            problems <- c(
              problems,
              sprintf(
                "L%d %s lists missing modified path: %s",
                line_no,
                session_id,
                target
              )
            )
          }
        }
      } else {
        problems <- c(
          problems,
          sprintf(
            "L%d %s linked log lacks ## Files Modified",
            line_no,
            session_id
          )
        )
      }
      if (
        any(grepl(
          "\\.log-scribe-[^[:space:]]*(auth|authentication|failed)",
          log_lines,
          ignore.case = TRUE,
          perl = TRUE
        ))
      ) {
        problems <- c(
          problems,
          sprintf(
            "L%d %s includes an authentication-failure artifact",
            line_no,
            session_id
          )
        )
      }
    }
    problems
  }

  fixture_root <- tempfile("espi-provenance-fixture-")
  fixture_log_dir <- file.path(fixture_root, ".living", "log")
  dir.create(fixture_log_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("known answer", file.path(fixture_root, "fixture.txt"))
  writeLines(
    c(
      "---",
      "session_id: fixture-valid",
      "ended: 2026-07-14T00:00:00-0500",
      "duration_minutes: 2",
      "files_changed: 1",
      "---",
      "",
      "## Session Summary",
      "valid",
      "",
      "## Key Outputs",
      "fixture",
      "",
      "## Status",
      "complete",
      "",
      "## Files Modified",
      "- fixture.txt"
    ),
    file.path(fixture_log_dir, "valid.md")
  )
  writeLines(
    c(
      "---",
      "session_id: fixture-malformed",
      "ended: NA",
      "---",
      "",
      "## Session Log"
    ),
    file.path(fixture_log_dir, "malformed.md")
  )
  fixture_registry <- file.path(fixture_log_dir, "LOG_REGISTRY.md")
  writeLines(
    c(
      "# Fixture",
      "| Date | Session ID | Project | Branch | Duration | Files Changed | Summary | Key Outputs | Status | Tags | Log |",
      "|------|------------|---------|--------|----------|---------------|---------|-------------|--------|------|-----|",
      "| 2026-07-14 | fixture-valid | espi | main | 2m | 1 | valid summary | fixture output | complete | fixture | [log](valid.md) |",
      "| 2026-07-14 | fixture-malformed | espi | main | 2m | 1 | malformed summary | fixture output | complete | fixture | [log](malformed.md) |"
    ),
    fixture_registry
  )
  fixture_problems <- validate_registry(
    fixture_registry,
    c("fixture-valid", "fixture-malformed"),
    c("fixture-valid", "fixture-malformed")
  )
  if (!any(grepl("fixture-malformed", fixture_problems, fixed = TRUE))) {
    return(fail(
      slug,
      "Known-answer malformed provenance fixture did not fail validation."
    ))
  }
  if (!file.exists(registry)) {
    return(fail(slug, ".living/log/LOG_REGISTRY.md is missing."))
  }
  registry_lines <- read_text(registry)
  table_lines <- grep("^\\|", registry_lines)
  header <- parse_markdown_table_row(registry_lines[[table_lines[[1L]]]])
  session_idx <- match("Session ID", header)
  status_idx <- match("Status", header)
  enforced_ids <- character()
  if (!is.na(session_idx) && !is.na(status_idx)) {
    rows <- lapply(table_lines[-REGISTRY_HEADER_ROW_INDEX], function(i) {
      fields <- parse_markdown_table_row(registry_lines[[i]])
      if (length(fields) != length(header)) {
        return(NULL)
      }
      stats::setNames(fields, header)
    })
    rows <- Filter(Negate(is.null), rows)
    ids <- vapply(rows, `[[`, character(1), "Session ID")
    statuses <- vapply(rows, `[[`, character(1), "Status")
    # Enforce the full semantic contract prospectively from the session that
    # introduced it; earlier hook-generated rows remain historical records.
    enforcement_start <- "2026-07-14-005"
    enforced_ids <- ids[
      grepl("^\\d{4}-\\d{2}-\\d{2}-\\d{3}$", ids) &
        ids >= enforcement_start &
        tolower(trimws(statuses)) == "complete"
    ]
  }
  scoped_ids <- unique(c("2026-07-05-007", "2026-07-05-008", enforced_ids))
  problems <- validate_registry(registry, scoped_ids, enforced_ids)
  if (length(problems) > 0L) {
    return(fail(slug, compact_problem_list(problems)))
  }
  pass(
    slug,
    sprintf(
      "Validated %d complete review/new registry row(s), linked logs, frontmatter consistency, required sections, modified paths, and authentication-artifact exclusion.",
      length(scoped_ids)
    )
  )
}


# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_metadata_contract <- function(root) {
  # Scientific boundary: Mouse, Condition, and derived sample_id define the
  # pseudobulk sample identity; missing or drifting metadata changes the biology.
  slug <- "metadata-contract"
  labels_path <- file.path(root, "analysis_labels.yml")
  preprocess_path <- file.path(root, "scripts", "03-preprocess.R")
  if (!file.exists(labels_path)) {
    return(fail(slug, "analysis_labels.yml is missing."))
  }
  if (!file.exists(preprocess_path)) {
    return(fail(slug, "scripts/03-preprocess.R is missing."))
  }

  yml <- read_text(labels_path)
  required <- extract_yaml_list(yml, "required_columns")
  required_ok <- all(c("Mouse", "Condition") %in% required)
  derived_ok <- any(grepl("^[[:space:]]*sample_id:", yml))

  lines <- read_text(preprocess_path)
  code <- squash(strip_inline_comment(lines))
  sample_id_ok <- grepl("sample_id[[:space:]]*<-", code) &&
    grepl("Mouse", code) &&
    grepl("Condition", code)
  has_validator_call <- grepl(
    "validate_required_metadata\\s*\\(.*Mouse.*Condition",
    code,
    perl = TRUE
  )
  validator_code <- ""
  validator_files <- list.files(
    file.path(root, "R"),
    pattern = "\\.[Rr]$",
    full.names = TRUE
  )
  for (validator_file in validator_files[file.exists(validator_files)]) {
    validator_code <- paste(validator_code, squash(read_text(validator_file)))
  }
  validator_checks_columns <- has_validator_call &&
    grepl("validate_required_metadata", validator_code, fixed = TRUE) &&
    grepl(
      "setdiff\\s*\\(\\s*columns\\s*,\\s*names\\s*\\(",
      validator_code,
      perl = TRUE
    )
  validator_checks_values <- has_validator_call &&
    grepl("validate_required_metadata", validator_code, fixed = TRUE) &&
    grepl("is\\.na\\s*\\(", validator_code, perl = TRUE) &&
    grepl("trimws\\s*\\(", validator_code, perl = TRUE)
  has_inline_column_guard <- grepl(
    "Mouse.*Condition.*(colnames|names|%in%|setdiff)|required_cols",
    code,
    ignore.case = TRUE,
    perl = TRUE
  )
  has_inline_missing_guard <- grepl(
    "(Mouse|Condition).*(is\\.na|anyNA|complete\\.cases|nzchar|trimws)",
    code,
    ignore.case = TRUE,
    perl = TRUE
  )
  has_column_guard <- validator_checks_columns || has_inline_column_guard
  has_missing_guard <- validator_checks_values || has_inline_missing_guard

  problems <- character()
  if (!required_ok) {
    problems <- c(
      problems,
      "analysis_labels.yml does not declare required_columns Mouse and Condition"
    )
  }
  if (!derived_ok) {
    problems <- c(
      problems,
      "analysis_labels.yml does not declare derived sample_id"
    )
  }
  if (!sample_id_ok) {
    problems <- c(
      problems,
      "preprocess-sobj.R does not derive sample_id from Mouse and Condition"
    )
  }
  if (!has_column_guard) {
    problems <- c(
      problems,
      "preprocess-sobj.R lacks a static guard that required metadata columns exist before sample_id construction"
    )
  }
  if (!has_missing_guard) {
    problems <- c(
      problems,
      "preprocess-sobj.R lacks a static guard for missing/blank Mouse or Condition values before sample_id construction"
    )
  }

  if (length(problems) > 0) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }
  pass(
    slug,
    "analysis_labels.yml declares Mouse, Condition, and derived sample_id; preprocess-sobj.R validates and derives them before downstream use."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
condition_in_blind_calls <- function(lines) {
  starts <- grep(
    "\\b(FindVariableFeatures|run_log1p_pca|run_pflog_pca|RunUMAP|FindClusters|FindNeighbors)\\s*\\(",
    lines,
    perl = TRUE
  )
  bad <- integer()
  for (start in starts) {
    depth <- 0L
    for (i in seq.int(start, length(lines))) {
      chars <- strsplit(lines[[i]], "", fixed = TRUE)[[1]]
      depth <- depth + sum(chars == "(") - sum(chars == ")")
      if (grepl("Condition", lines[[i]], fixed = TRUE)) {
        bad <- c(bad, i)
      }
      if (depth <= 0L) {
        break
      }
    }
  }
  unique(bad)
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_label_firewall <- function(root) {
  slug <- "label-permutation"
  labels_path <- file.path(root, "analysis_labels.yml")
  if (!file.exists(labels_path)) {
    return(fail(
      slug,
      "analysis_labels.yml is missing, so the label firewall has no declared label column."
    ))
  }
  yml <- read_text(labels_path)
  if (!identical(extract_yaml_scalar(yml, "treatment"), "Condition")) {
    return(fail(
      slug,
      "analysis_labels.yml must declare labels.treatment as Condition."
    ))
  }
  paths <- file.path(root, "scripts", c("03-preprocess.R", "04-cluster.R"))
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    return(fail(
      slug,
      paste(
        "Missing blind-stage script(s):",
        paste(basename(missing), collapse = ", ")
      )
    ))
  }
  bad_refs <- character()
  for (path in paths) {
    lines <- strip_inline_comment(read_text(path))
    blind_bad <- condition_in_blind_calls(lines)
    if (length(blind_bad) > 0L) {
      bad_refs <- c(
        bad_refs,
        paste0(basename(path), " ", line_refs(lines, blind_bad))
      )
    }
    condition_lines <- which(grepl("Condition", lines, fixed = TRUE))
    allowed <- grepl(
      "sample_id|metadata|required|colnames|missing|anyNA|complete\\.cases|is\\.na|nzchar|label|treatment|paste0|gsub|permute|sample_table|\"Condition\"|sobj\\$Condition",
      lines,
      ignore.case = TRUE,
      perl = TRUE
    )
    general_bad <- setdiff(
      condition_lines[!allowed[condition_lines]],
      blind_bad
    )
    if (length(general_bad) > 0L) {
      bad_refs <- c(
        bad_refs,
        paste0(basename(path), " ", line_refs(lines, general_bad))
      )
    }
  }
  if (length(bad_refs) > 0L) {
    return(fail(
      slug,
      paste(
        "Condition appears to influence blind-stage code outside metadata/sample_id boundaries:",
        paste(bad_refs, collapse = " | ")
      )
    ))
  }

  path_env <- new.env(parent = globalenv())
  path_error <- tryCatch(
    {
      sys.source(file.path(root, "R", "paths.R"), envir = path_env)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (
    !is.null(path_error) ||
      !exists("INPUT_OBJECT_DIR", envir = path_env, inherits = FALSE)
  ) {
    return(skip(
      slug,
      "Static firewall passed, but production path resolution is unavailable for a scratch permutation run."
    ))
  }
  input_path <- file.path(
    get("INPUT_OBJECT_DIR", envir = path_env),
    "sobj_qc_filtered.rds"
  )
  if (!file.exists(input_path)) {
    return(skip(
      slug,
      "Static firewall passed, but INPUT_OBJECT_DIR/sobj_qc_filtered.rds is unavailable for a scratch permutation run."
    ))
  }
  rscript <- Sys.which("Rscript")
  script <- file.path(root, "scripts", "03-preprocess.R")
  if (identical(unname(rscript), "")) {
    return(fail(
      slug,
      "Rscript is unavailable, so label permutation cannot be executed."
    ))
  }
  scratch <- tempfile("espi-label-permutation-")
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)
  run_dirs <- file.path(scratch, c("baseline", "permuted"))
  for (run_dir in run_dirs) {
    dir.create(
      file.path(run_dir, "seurat_objects", "input"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (
      !file.copy(
        input_path,
        file.path(run_dir, "seurat_objects", "input", "sobj_qc_filtered.rds")
      )
    ) {
      return(skip(
        slug,
        "QC-filtered input could not be copied into both scratch workspaces."
      ))
    }
  }
  # ANALYSIS_OK[file-freshness-tripwire]: mtime snapshots are the explicit side-effect signal; the before/after comparison below detects unexpected production writes.
  snapshots <- function(paths) {
    out <- character()
    for (path in paths) {
      if (!dir.exists(path)) {
        next
      }
      files <- list.files(
        path,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
      )
      if (length(files) == 0L) {
        next
      }
      info <- file.info(files)
      keys <- normalizePath(files, winslash = "/", mustWork = FALSE)
      out <- c(out, stats::setNames(paste(info$size, info$mtime), keys))
    }
    out
  }
  production_dirs <- c(
    get("CURRENT_OBJECT_DIR", envir = path_env),
    get("FIGURE_DIR", envir = path_env)
  )
  before <- snapshots(production_dirs)
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  run_once <- function(run_dir, checkpoint_path, seed = NULL) {
    names_to_set <- c(
      "DATA_ROOT_DIR",
      "ESPI_TRIPWIRE_MODE",
      "CHECKPOINT_LOG",
      "STOP_AFTER_CHECKPOINT"
    )
    old <- Sys.getenv(
      c(names_to_set, "ESPI_PERMUTE_CONDITION_SEED"),
      unset = NA_character_
    )
    on.exit(
      {
        for (i in seq_along(names_to_set)) {
          if (is.na(old[[i]])) {
            Sys.unsetenv(names_to_set[[i]])
          } else {
            Sys.setenv(structure(old[[i]], names = names_to_set[[i]]))
          }
        }
        if (is.na(old[[length(names_to_set) + 1L]])) {
          Sys.unsetenv("ESPI_PERMUTE_CONDITION_SEED")
        } else {
          Sys.setenv(
            ESPI_PERMUTE_CONDITION_SEED = old[[length(names_to_set) + 1L]]
          )
        }
      },
      add = TRUE
    )
    Sys.setenv(
      DATA_ROOT_DIR = run_dir,
      ESPI_TRIPWIRE_MODE = "true",
      CHECKPOINT_LOG = checkpoint_path,
      STOP_AFTER_CHECKPOINT = "blind_qc_complete"
    )
    if (is.null(seed)) {
      Sys.unsetenv("ESPI_PERMUTE_CONDITION_SEED")
    } else {
      Sys.setenv(ESPI_PERMUTE_CONDITION_SEED = as.character(seed))
    }
    tryCatch(
      system2(
        rscript,
        c(
          script,
          "--input",
          file.path(run_dir, "seurat_objects", "input", "sobj_qc_filtered.rds"),
          "--normalization",
          "log1p"
        ),
        stdout = TRUE,
        stderr = TRUE
      ),
      warning = function(w) structure(conditionMessage(w), status = 1L),
      error = function(e) structure(conditionMessage(e), status = 1L)
    )
  }
  baseline_log <- file.path(scratch, "baseline.tsv")
  permuted_log <- file.path(scratch, "permuted.tsv")
  baseline_output <- run_once(run_dirs[[1L]], baseline_log)
  # ANALYSIS_OK[positional-fixture-index]: run_dirs is constructed as baseline then permuted immediately above; the explicit index is checked by paired output comparisons below.
  permuted_output <- run_once(
    run_dirs[[PERMUTED_RUN_INDEX]],
    permuted_log,
    seed = 20260714L
  )
  after <- snapshots(production_dirs)
  if (!identical(before, after)) {
    return(fail(
      slug,
      "Label permutation changed files beneath production CURRENT_OBJECT_DIR or FIGURE_DIR."
    ))
  }
  # ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
  status_of <- function(output) {
    status <- attr(output, "status")
    if (is.null(status)) 0L else as.integer(status)
  }
  if (status_of(baseline_output) != 0L || status_of(permuted_output) != 0L) {
    return(fail(
      slug,
      sprintf(
        "Scratch preprocessing failed (baseline exit=%d, permuted exit=%d).",
        status_of(baseline_output),
        status_of(permuted_output)
      )
    ))
  }
  # ANALYSIS_OK[optional-input]: missing or malformed fingerprint logs are intentional fault-injection inputs; callers assert controlled failure rather than treating them as valid output.
  parse_fingerprint <- function(path) {
    if (!file.exists(path)) {
      return(NULL)
    }
    log <- tryCatch(
      utils::read.delim(
        path,
        sep = "\t",
        quote = "\"",
        stringsAsFactors = FALSE
      ),
      error = function(e) {
        message(sprintf("Fingerprint log read failed: %s", conditionMessage(e)))
        data.frame()
      }
    )
    if (
      is.null(log) || !"checkpoint" %in% names(log) || !"fields" %in% names(log)
    ) {
      return(NULL)
    }
    rows <- log[log$checkpoint == "blind_qc_complete", , drop = FALSE]
    if (nrow(rows) != 1L) {
      return(NULL)
    }
    fields <- strsplit(as.character(rows$fields[[1L]]), ";", fixed = TRUE)[[1L]]
    values <- strsplit(fields, "=", fixed = TRUE)
    parsed <- stats::setNames(
      vapply(values, function(x) paste(x[-1L], collapse = "="), character(1)),
      vapply(values, `[[`, character(1), 1L)
    )
    parsed
  }
  baseline <- parse_fingerprint(baseline_log)
  permuted <- parse_fingerprint(permuted_log)
  required <- c(
    "fingerprint_algorithm",
    "hvg_fingerprint",
    "pca_sdev_fingerprint"
  )
  if (
    is.null(baseline) ||
      is.null(permuted) ||
      any(!required %in% names(baseline)) ||
      any(!required %in% names(permuted))
  ) {
    return(fail(
      slug,
      "Both scratch runs must emit one blind_qc_complete checkpoint with exact fingerprint fields."
    ))
  }
  if (
    !identical(baseline[["fingerprint_algorithm"]], "exact-v1") ||
      !identical(permuted[["fingerprint_algorithm"]], "exact-v1")
  ) {
    return(fail(
      slug,
      "blind_qc_complete must declare fingerprint_algorithm=exact-v1."
    ))
  }
  if (
    !identical(baseline[["hvg_fingerprint"]], permuted[["hvg_fingerprint"]]) ||
      !identical(
        baseline[["pca_sdev_fingerprint"]],
        permuted[["pca_sdev_fingerprint"]]
      )
  ) {
    return(fail(
      slug,
      "Condition permutation changed the declared blind HVG or PCA fingerprints."
    ))
  }
  pass(
    slug,
    "Static firewall passed and original/permuted scratch runs produced identical exact-v1 HVG and PCA fingerprints without production writes."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
tripwire_toy_contrast_direction <- function(root) {
  # Scientific boundary: future DE signs must mean target/reference, not the
  # accidental reverse comparison.
  slug <- "toy-contrast-direction"
  labels_path <- file.path(root, "analysis_labels.yml")
  if (!file.exists(labels_path)) {
    return(fail(slug, "analysis_labels.yml is missing."))
  }
  yml <- read_text(labels_path)
  target <- extract_yaml_scalar(yml, "target")
  reference <- extract_yaml_scalar(yml, "reference")
  sign <- extract_yaml_scalar(yml, "sign_convention")
  unit <- extract_yaml_scalar(yml, "statistical_unit")

  problems <- character()
  if (!identical(target, "p27CKO +EStim")) {
    problems <- c(problems, "contrast target must be p27CKO +EStim")
  }
  if (!identical(reference, "p27CKO")) {
    problems <- c(problems, "contrast reference must be p27CKO")
  }
  if (
    is.na(sign) ||
      !grepl("higher in p27CKO \\+EStim than p27CKO", sign, fixed = FALSE)
  ) {
    problems <- c(
      problems,
      "sign convention must define positive log fold-change as higher in p27CKO +EStim than p27CKO"
    )
  }
  if (
    is.na(unit) ||
      !grepl("Mouse.*Condition.*pseudobulk", unit, ignore.case = TRUE)
  ) {
    problems <- c(
      problems,
      "statistical unit must be Mouse x Condition pseudobulk sample"
    )
  }
  if (length(problems) > 0) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }

  code_paths <- c(
    list.files(
      file.path(root, "scripts"),
      pattern = "\\.[Rr]$",
      full.names = TRUE
    ),
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
  )
  existing_code_paths <- code_paths[file.exists(code_paths)]
  de_entry_pattern <- paste(
    "\\b(FindMarkers|FindAllMarkers|DESeq2|edgeR|limma)\\b",
    "\\bdifferential[ -]?expression\\b",
    "\\b(logFC|log2FoldChange)\\b",
    sep = "|"
  )
  has_de_entry <- FALSE
  if (length(existing_code_paths) > 0) {
    for (path in existing_code_paths) {
      txt <- squash(read_text(path))
      if (grepl(de_entry_pattern, txt, ignore.case = TRUE, perl = TRUE)) {
        has_de_entry <- TRUE
        break
      }
    }
  }

  if (!has_de_entry) {
    return(skip(
      slug,
      "Contrast direction is encoded in analysis_labels.yml; toy known-answer DE run is skipped because no differential-expression entry point exists yet."
    ))
  }

  de_path <- file.path(root, "scripts", "12-run-mg-de.R")
  if (file.exists(de_path)) {
    de_txt <- squash(read_text(de_path))
    required_patterns <- c(
      'CONTROL_LEVEL <- "control"',
      'ESTIM_LEVEL <- "estim"',
      'CONTRAST_DIRECTION <- "estim_vs_control"',
      'levels = c(CONTROL_LEVEL, ESTIM_LEVEL)',
      'contrast = c("condition", ESTIM_LEVEL, CONTROL_LEVEL)'
    )
    missing_patterns <- required_patterns[
      !vapply(
        required_patterns,
        function(pattern) grepl(pattern, de_txt, fixed = TRUE),
        logical(1)
      )
    ]
    if (length(missing_patterns) == 0) {
      return(pass(
        slug,
        "run-mg-selected-de.R encodes estim/control factor order and DESeq2 contrast direction."
      ))
    }
    return(fail(
      slug,
      paste(
        "run-mg-selected-de.R is missing contrast-direction guard(s):",
        paste(missing_patterns, collapse = "; ")
      )
    ))
  }

  fail(
    slug,
    "A differential-expression entry point appears to exist, but this runner has no toy known-answer execution hook for it yet."
  )
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
print_results <- function(results) {
  status_width <- max(
    nchar("status"),
    nchar(vapply(results, function(x) x$status, character(1)))
  )
  slug_width <- max(
    nchar("slug"),
    nchar(vapply(results, function(x) x$slug, character(1)))
  )
  header <- sprintf(
    "%-*s  %-*s  %s",
    status_width,
    "status",
    slug_width,
    "slug",
    "message"
  )
  sep <- sprintf(
    "%-*s  %-*s  %s",
    status_width,
    paste(rep("-", status_width), collapse = ""),
    slug_width,
    paste(rep("-", slug_width), collapse = ""),
    paste(rep("-", 7), collapse = "")
  )
  writeLines(header)
  writeLines(sep)
  for (res in results) {
    writeLines(sprintf(
      "%-*s  %-*s  %s",
      status_width,
      res$status,
      slug_width,
      res$slug,
      res$message
    ))
  }
}

# ANALYSIS_OK[script-entrypoint]: internal helper is dispatched by main() in this executable; R026 does not model same-file entrypoints, and main() invokes the check.
main <- function() {
  root <- find_repo_root()
  if (is.na(root)) {
    results <- list(fail(
      "runner-root",
      "Run from the ESPI repo root or invoke tools/run-tripwires.R from inside the repo."
    ))
    print_results(results)
    quit(status = 1L, save = "no")
  }
  setwd(root)

  checks <- list(
    tripwire_cluster_wrapper_contract,
    tripwire_cli_value_boundaries,
    tripwire_pipeline_dry_run_contract,
    tripwire_just_public_interface,
    tripwire_branch_artifact_collision,
    tripwire_report_values_freshness,
    tripwire_missing_counts_file,
    tripwire_missing_metadata_sample,
    tripwire_heatmap_missing_input,
    tripwire_p27_rng_state_preservation,
    tripwire_mycelium_provenance_semantics,
    tripwire_metadata_contract,
    tripwire_label_firewall,
    tripwire_toy_contrast_direction
  )

  results <- lapply(checks, function(fun) {
    tryCatch(fun(root), error = function(e) {
      fail("runner-error", conditionMessage(e))
    })
  })
  print_results(results)
  has_fail <- any(vapply(
    results,
    function(x) identical(x$status, "FAIL"),
    logical(1)
  ))
  quit(status = if (has_fail) 1L else 0L, save = "no")
}

main()
