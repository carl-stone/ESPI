#' Three-color analysis palette for depleted/neutral/enriched (low/mid/high) encodings.
#' @export
palette_analysis_three <- c(low = "#2166ac", mid = "grey75", high = "#e31a8c")

#' Two-color low/high analysis palette; subset of `palette_analysis_three`.
#' @export
palette_dotplot_pair <- unname(palette_analysis_three[c("low", "high")])
