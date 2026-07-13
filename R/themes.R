#' Three-color analysis palette for depleted/neutral/enriched (low/mid/high) encodings.
#' @export
palette_analysis_three <- c(low = "#2166ac", mid = "grey75", high = "#e31a8c")

#' Two-color low/high analysis palette; subset of `palette_analysis_three`.
#' @export
palette_dotplot_pair <- unname(palette_analysis_three[c("low", "high")])

#' Theme for all figures
#' @export
theme_stone <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      axis.text = ggplot2::element_text(face = "bold", color = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
}
