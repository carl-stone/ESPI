# conditions.R — Condition metadata labels and display strings.

#' Metadata column containing condition labels.
#' @export
CONDITION_COL <- "Condition"
#' E-Stim condition label (metadata value).
#' @export
ESTIM_LABEL <- "p27CKO +EStim"
#' Control condition label (metadata value).
#' @export
CTRL_LABEL <- "p27CKO"

#' Display label for the E-Stim condition.
#' @export
ESTIM_DISPLAY_LABEL <- "p27CKO + E-Stim"
#' Display label for the control condition.
#' @export
CTRL_DISPLAY_LABEL <- "p27CKO"
#' Parenthetical contrast label for plot axes/titles.
#' @export
CONTRAST_DISPLAY_LABEL <- sprintf(
  "(%s vs. %s)",
  ESTIM_DISPLAY_LABEL,
  CTRL_DISPLAY_LABEL
)
