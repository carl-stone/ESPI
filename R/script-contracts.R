#' Order labels deterministically.
#'
#' Numeric labels sort numerically; all other labels sort lexically.
#'
#' @param values Label values.
#'
#' @return A character vector of unique ordered labels.
#' @export
cluster_levels_for_labels <- function(values) {
  labels <- unique(as.character(values))
  if (all(grepl("^-?[0-9]+$", labels))) {
    return(as.character(sort(as.integer(labels))))
  }
  sort(labels, method = "radix")
}

#' Require a non-empty character scalar.
#'
#' @param x Value to validate.
#' @param name Argument name for error messages.
#'
#' @return `x`, invisibly.
#' @export
assert_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(name, " must be a non-empty character scalar.", call. = FALSE)
  }
  invisible(x)
}
