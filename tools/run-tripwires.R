#!/usr/bin/env Rscript

# ESPI tripwire runner.
# All checks use base R and avoid running the full analysis or rendering Quarto.

make_result <- function(status, slug, message) {
  list(status = status, slug = slug, message = message)
}

pass <- function(slug, message) make_result("PASS", slug, message)
fail <- function(slug, message) make_result("FAIL", slug, message)
skip <- function(slug, message) make_result("SKIP", slug, message)

script_path <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  sub("^--file=", "", file_arg[[1]])
}

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

compact_problem_list <- function(problems, max_n = 8L) {
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

parse_markdown_table_row <- function(line) {
  line <- sub("^\\|", "", line)
  line <- sub("\\|[[:space:]]*$", "", line)
  fields <- trimws(strsplit(line, "\\|", fixed = FALSE)[[1]])
  if (length(fields) == 0L) {
    return(character())
  }
  fields
}

extract_markdown_link_target <- function(text) {
  match <- regexec("\\[[^]]+\\]\\(([^)]+)\\)", text, perl = TRUE)
  hit <- regmatches(text, match)[[1]]
  if (length(hit) < 2L) {
    return(NA_character_)
  }
  hit[[2L]]
}

frontmatter_scalar <- function(lines, key) {
  if (length(lines) < 3L || !identical(lines[[1L]], "---")) {
    return(NA_character_)
  }
  end <- which(lines[-1L] == "---")
  if (length(end) == 0L) {
    return(NA_character_)
  }
  frontmatter <- lines[seq.int(2L, end[[1L]])]
  extract_yaml_scalar(frontmatter, key)
}

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
    Mouse = rep(c("M1", "M2", "M3"), each = 4L),
    Condition = rep(c("control", "estim", "control"), each = 4L),
    cluster = rep(c("1", "1", "2", "2"), times = 3L),
    row.names = cells,
    stringsAsFactors = FALSE
  )
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

restore_random_seed <- function(seed_state, existed) {
  if (existed) {
    assign(".Random.seed", seed_state, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

analysis_table_annotation_dir <- function(root) {
  paths_env <- new.env(parent = .GlobalEnv)
  sys.source(file.path(root, "R", "paths.R"), envir = paths_env)
  file.path(get("TABLE_DIR", envir = paths_env, inherits = FALSE), "annotation")
}


MAX_LINE_REFS <- 6L
YAML_LIST_CAPTURE_LENGTH <- 2L
YAML_LIST_CAPTURE_VALUE <- 2L


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

strip_inline_comment <- function(line) {
  # Good enough for the current scripts: ignore comment-only prose without trying
  # to parse quoted '#'. The tripwires below search code-like tokens.
  sub("^[[:space:]]*#.*$", "", line)
}

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

tripwire_mycelium_provenance_semantics <- function(root) {
  # Provenance boundary: review F5 is scoped to the heatmap session records.
  # Broader historical Mycelium registry backfill is a separate task, so older
  # blank/file-list-only rows outside this review must not fail this tripwire.
  slug <- "mycelium-provenance-semantics"
  registry <- file.path(root, ".living", "log", "LOG_REGISTRY.md")
  scoped_session_ids <- c("2026-07-05-007", "2026-07-05-008")
  if (!file.exists(registry)) {
    return(fail(slug, ".living/log/LOG_REGISTRY.md is missing."))
  }

  lines <- read_text(registry)
  table_lines <- grep("^\\|", lines)
  if (length(table_lines) < 3L) {
    return(fail(slug, "LOG_REGISTRY.md does not contain a markdown table."))
  }
  header <- parse_markdown_table_row(lines[[table_lines[[1L]]]])
  required <- c("Session ID", "Summary", "Key Outputs", "Status", "Tags", "Log")
  missing_columns <- setdiff(required, header)
  if (length(missing_columns) > 0L) {
    return(fail(
      slug,
      paste(
        "LOG_REGISTRY.md is missing required column(s):",
        paste(missing_columns, collapse = ", ")
      )
    ))
  }

  registry_rows <- list()
  problems <- character()
  for (line_no in table_lines[-seq_len(2L)]) {
    fields <- parse_markdown_table_row(lines[[line_no]])
    if (length(fields) != length(header)) {
      session_idx <- match("Session ID", header)
      session_id <- if (length(fields) >= session_idx) {
        fields[[session_idx]]
      } else {
        ""
      }
      if (length(session_id) == 1L && session_id %in% scoped_session_ids) {
        problems <- c(
          problems,
          sprintf(
            "L%d %s malformed registry row has %d fields",
            line_no,
            session_id,
            length(fields)
          )
        )
      }
      next
    }
    row <- stats::setNames(fields, header)
    if (row[["Session ID"]] %in% scoped_session_ids) {
      registry_rows[[row[["Session ID"]]]] <- list(row = row, line_no = line_no)
    }
  }

  missing_rows <- setdiff(scoped_session_ids, names(registry_rows))
  if (length(missing_rows) > 0L) {
    problems <- c(
      problems,
      sprintf(
        "LOG_REGISTRY.md is missing review-scoped row(s): %s",
        paste(missing_rows, collapse = ", ")
      )
    )
  }

  for (row_id in intersect(scoped_session_ids, names(registry_rows))) {
    entry <- registry_rows[[row_id]]
    row <- entry$row
    line_no <- entry$line_no

    if (!identical(tolower(trimws(row[["Status"]])), "complete")) {
      problems <- c(
        problems,
        sprintf("L%d %s status is not complete", line_no, row_id)
      )
    }
    for (field in c("Summary", "Key Outputs", "Tags")) {
      value <- trimws(row[[field]])
      if (!nzchar(value) || identical(value, "—")) {
        problems <- c(
          problems,
          sprintf("L%d %s has empty %s", line_no, row_id, field)
        )
      }
    }
    if (is_file_list_only_summary(row[["Summary"]])) {
      problems <- c(
        problems,
        sprintf("L%d %s summary is file-list-only", line_no, row_id)
      )
    }

    log_target <- extract_markdown_link_target(row[["Log"]])
    if (is.na(log_target) || !nzchar(log_target)) {
      problems <- c(
        problems,
        sprintf("L%d %s lacks a linked log", line_no, row_id)
      )
      next
    }
    log_path <- file.path(dirname(registry), log_target)
    if (!file.exists(log_path)) {
      problems <- c(
        problems,
        sprintf("L%d %s linked log is missing: %s", line_no, row_id, log_target)
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
            row_id,
            field
          )
        )
      }
    }
  }

  if (length(problems) > 0L) {
    return(fail(slug, compact_problem_list(problems)))
  }

  pass(
    slug,
    "Review-scoped LOG_REGISTRY rows 2026-07-05-007 and 2026-07-05-008 have semantic Summary/Key Outputs/Tags fields, and linked logs carry ended/duration/files_changed frontmatter."
  )
}


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

tripwire_label_firewall <- function(root) {
  # Scientific boundary: Condition is an interpretation label, not an input to
  # blind HVG selection, PCA, UMAP, or clustering.
  slug <- "label-permutation"
  labels_path <- file.path(root, "analysis_labels.yml")
  if (!file.exists(labels_path)) {
    return(fail(
      slug,
      "analysis_labels.yml is missing, so the label firewall has no declared label column."
    ))
  }
  yml <- read_text(labels_path)
  treatment <- extract_yaml_scalar(yml, "treatment")
  if (!identical(treatment, "Condition")) {
    return(fail(
      slug,
      "analysis_labels.yml must declare labels.treatment as Condition."
    ))
  }

  paths <- file.path(root, "scripts", c("03-preprocess.R", "04-cluster.R"))
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
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
    if (length(blind_bad) > 0) {
      bad_refs <- c(
        bad_refs,
        paste0(basename(path), " ", line_refs(lines, blind_bad))
      )
    }

    condition_lines <- which(grepl("Condition", lines, fixed = TRUE))
    allowed <- grepl(
      "sample_id|metadata|required|colnames|missing|anyNA|complete\\.cases|is\\.na|nzchar|label|treatment|paste0|gsub|sobj\\$Condition",
      lines,
      ignore.case = TRUE,
      perl = TRUE
    )
    general_bad <- setdiff(
      condition_lines[!allowed[condition_lines]],
      blind_bad
    )
    if (length(general_bad) > 0) {
      bad_refs <- c(
        bad_refs,
        paste0(basename(path), " ", line_refs(lines, general_bad))
      )
    }
  }

  if (length(bad_refs) > 0) {
    return(fail(
      slug,
      paste(
        "Condition appears to influence blind-stage code outside metadata/sample_id boundaries:",
        paste(bad_refs, collapse = " | ")
      )
    ))
  }

  skip(
    slug,
    "Static firewall passed: Condition is confined to metadata/sample_id boundaries in early scripts. Full label permutation is skipped because the project has no scratch-output hook for safe blind HVG/PCA/UMAP reruns."
  )
}

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
    tripwire_branch_artifact_collision,
    tripwire_report_values_freshness,
    tripwire_missing_counts_file,
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
