#' Emit a tripwire checkpoint.
#'
#' When `CHECKPOINT_LOG` is set, append one tab-separated row describing a named
#' scientific boundary. When `STOP_AFTER_CHECKPOINT` equals `checkpoint`, stop
#' the current R process after the row is written.
#'
#' @param checkpoint Character scalar checkpoint name.
#' @param ... Named scalar fields to record with the checkpoint.
#'
#' @return `NULL`, invisibly.
#' @export
# ANALYSIS_OK[R026]: package export is loaded by devtools::load_all and invoked by executable analysis scripts.
emit_tripwire_checkpoint <- function(checkpoint, ...) {
  stopifnot(
    length(checkpoint) == 1,
    is.character(checkpoint),
    !is.na(checkpoint)
  )

  fields <- list(...)
  log_path <- Sys.getenv("CHECKPOINT_LOG", unset = "")
  if (!identical(log_path, "")) {
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    log_exists <- file.exists(log_path)
    row <- data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
      checkpoint = checkpoint,
      status = "passed",
      fields = encode_checkpoint_fields(fields),
      stringsAsFactors = FALSE
    )
    utils::write.table(
      row,
      file = log_path,
      sep = "\t",
      row.names = FALSE,
      col.names = !log_exists,
      append = log_exists,
      quote = TRUE
    )
  }

  stop_after <- Sys.getenv("STOP_AFTER_CHECKPOINT", unset = "")
  if (identical(stop_after, checkpoint)) {
    quit(save = "no", status = 0)
  }

  invisible(NULL)
}

#' Validate required metadata columns.
#'
#' @param meta Data frame-like metadata table.
#' @param columns Character vector of required metadata columns.
#' @param stage Character scalar checkpoint stage name.
#'
#' @return `NULL`, invisibly.
#' @export
validate_required_metadata <- function(
  meta,
  columns,
  # ANALYSIS_OK[smuggled-default]: intentional function signature default for metadata stage.
  stage = "metadata_complete"
) {
  stopifnot(is.character(columns), length(columns) > 0)

  missing_cols <- setdiff(columns, names(meta))
  if (length(missing_cols) > 0) {
    write_tripwire_drop_ledger(
      sample_ids = character(),
      stage = stage,
      reason = paste0(
        "missing required columns: ",
        paste(missing_cols, collapse = ", ")
      ),
      allowed = FALSE
    )
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  missing_value_cols <- columns[vapply(
    columns,
    function(column) {
      values <- meta[[column]]
      any(is.na(values) | trimws(as.character(values)) == "")
    },
    logical(1)
  )]
  if (length(missing_value_cols) > 0) {
    bad_rows <- Reduce(
      union,
      lapply(missing_value_cols, function(column) {
        values <- meta[[column]]
        which(is.na(values) | trimws(as.character(values)) == "")
      })
    )
    sample_ids <- rownames(meta)[bad_rows]
    if (is.null(sample_ids)) {
      sample_ids <- as.character(bad_rows)
    }
    write_tripwire_drop_ledger(
      sample_ids = sample_ids,
      stage = stage,
      reason = paste0(
        "missing values in: ",
        paste(missing_value_cols, collapse = ", ")
      ),
      allowed = FALSE
    )
    stop(
      "Missing metadata values in required columns: ",
      paste(missing_value_cols, collapse = ", "),
      call. = FALSE
    )
  }

  emit_tripwire_checkpoint(
    stage,
    required_columns = paste(columns, collapse = ","),
    n_rows = nrow(meta)
  )

  invisible(NULL)
}

#' Record dropped samples or cells for tripwire audits.
#'
#' When `DROP_LEDGER` is unset, this function is a no-op.
#'
#' @param sample_ids Character vector of sample/cell identifiers.
#' @param stage Character scalar pipeline stage.
#' @param reason Character scalar drop reason.
#' @param allowed Logical scalar indicating whether the drop is expected policy.
#'
#' @return `NULL`, invisibly.
#' @export
# ANALYSIS_OK[R026]: package export is loaded by devtools::load_all and invoked by executable analysis scripts.
write_tripwire_drop_ledger <- function(
  sample_ids,
  stage,
  reason,
  allowed = FALSE
) {
  stopifnot(length(stage) == 1, length(reason) == 1, length(allowed) == 1)

  ledger_path <- Sys.getenv("DROP_LEDGER", unset = "")
  if (identical(ledger_path, "")) {
    return(invisible(NULL))
  }

  dir.create(dirname(ledger_path), recursive = TRUE, showWarnings = FALSE)
  ledger_exists <- file.exists(ledger_path)
  if (length(sample_ids) == 0) {
    sample_ids <- NA_character_
  }
  rows <- data.frame(
    sample_id = as.character(sample_ids),
    stage = stage,
    reason = reason,
    allowed_by_policy = isTRUE(allowed),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    rows,
    file = ledger_path,
    sep = "\t",
    row.names = FALSE,
    col.names = !ledger_exists,
    append = ledger_exists,
    quote = TRUE
  )

  invisible(NULL)
}

# ANALYSIS_OK[R026]: package helper is loaded by devtools::load_all and called by same-file checkpoint encoding.
encode_checkpoint_fields <- function(fields) {
  if (length(fields) == 0) {
    return("")
  }
  values <- vapply(fields, encode_checkpoint_value, character(1))
  paste(paste0(names(values), "=", values), collapse = ";")
}

# ANALYSIS_OK[R026]: package helper is loaded by devtools::load_all and called by same-file checkpoint encoding.
encode_checkpoint_value <- function(value) {
  if (length(value) == 0 || is.null(value)) {
    return("")
  }
  value <- paste(as.character(value), collapse = ",")
  gsub("[\t\r\n;]", " ", value)
}
