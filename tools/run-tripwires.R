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
    file.path("scripts", "cluster-sobj.R"),
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

line_refs <- function(lines, idx, max_n = 6) {
  if (length(idx) == 0) {
    return("")
  }
  idx <- idx[seq_len(min(length(idx), max_n))]
  paste(sprintf("L%d: %s", idx, trimws(lines[idx])), collapse = "; ")
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
    if (length(hit) == 2) {
      out <- c(out, hit[[2]])
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
  # overwrite, share, or masquerade as the same clustering artifact.
  slug <- "branch-artifact-collision"
  path <- file.path(root, "scripts", "cluster-sobj.R")
  if (!file.exists(path)) {
    return(fail(slug, "scripts/cluster-sobj.R is missing."))
  }
  lines <- read_text(path)
  code <- squash(strip_inline_comment(lines))

  has_cc_tag <- grepl(
    "cc_tag[[:space:]]*<-.*filtered_cell_cycle",
    code,
    perl = TRUE
  )
  has_branch_tag <- grepl(
    "branch_tag[[:space:]]*<-[[:space:]]*sprintf\\([[:space:]]*\"%s_%s\"[[:space:]]*,[[:space:]]*norm[[:space:]]*,[[:space:]]*cc_tag",
    code,
    perl = TRUE
  )

  persistent_name_line <- grepl(
    "sprintf\\([[:space:]]*\"(cluster_|umap_|%s_dims|.*\\.rds)",
    lines,
    perl = TRUE
  ) |
    grepl(
      "(prefix|out_tag|reduction_name|color_by)[[:space:]]*=",
      lines,
      perl = TRUE
    )
  norm_only <- which(
    persistent_name_line &
      grepl("\\bnorm\\b", lines) &
      !grepl("\\b(branch_tag|cc_tag)\\b", lines)
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
      "missing branch_tag <- sprintf(\"%s_%s\", norm, cc_tag)"
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

  if (length(problems) > 0) {
    return(fail(slug, paste(problems, collapse = " | ")))
  }
  pass(
    slug,
    "cluster columns, UMAP reductions, clustree tags, and clustered RDS names use the normalization + cell-cycle branch tag contract."
  )
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
  missing_figures <- figure_paths[!file.exists(figure_paths)]

  problems <- character()
  if (length(missing_figures) > 0) {
    problems <- c(
      problems,
      paste(
        "missing QMD figure reference(s):",
        paste(
          file.path("notebook", figure_refs[!file.exists(figure_paths)]),
          collapse = ", "
        )
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
  script <- file.path(root, "scripts", "preprocess-sobj.R")
  if (!file.exists(script)) {
    return(fail(slug, "scripts/preprocess-sobj.R is missing."))
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

tripwire_metadata_contract <- function(root) {
  # Scientific boundary: Mouse, Condition, and derived sample_id define the
  # pseudobulk sample identity; missing or drifting metadata changes the biology.
  slug <- "metadata-contract"
  labels_path <- file.path(root, "analysis_labels.yml")
  preprocess_path <- file.path(root, "scripts", "preprocess-sobj.R")
  if (!file.exists(labels_path)) {
    return(fail(slug, "analysis_labels.yml is missing."))
  }
  if (!file.exists(preprocess_path)) {
    return(fail(slug, "scripts/preprocess-sobj.R is missing."))
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

  paths <- file.path(root, "scripts", c("preprocess-sobj.R", "cluster-sobj.R"))
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
  code_paths <- code_paths[file.exists(code_paths)]
  de_entry_pattern <- paste(
    "\\b(FindMarkers|FindAllMarkers|DESeq2|edgeR|limma)\\b",
    "\\bdifferential[ -]?expression\\b",
    "\\b(logFC|log2FoldChange)\\b",
    sep = "|"
  )
  has_de_entry <- FALSE
  if (length(code_paths) > 0) {
    for (path in code_paths) {
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
    tripwire_branch_artifact_collision,
    tripwire_report_values_freshness,
    tripwire_missing_counts_file,
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
