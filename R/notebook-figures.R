#' Link a PNG into the analysis notebook figure directory.
#'
#' Replaces an existing file or symlink of the same name.
#'
#' @param png_path Path to the PNG artifact.
#'
#' @return The notebook symlink path.
#' @export
link_notebook_png <- function(png_path) {
  notebook_figure_dir <- here::here("notebook", "figures")
  dir.create(notebook_figure_dir, recursive = TRUE, showWarnings = FALSE)
  notebook_png_path <- file.path(notebook_figure_dir, basename(png_path))
  if (
    file.exists(notebook_png_path) || nzchar(Sys.readlink(notebook_png_path))
  ) {
    unlink(notebook_png_path)
  }
  if (!isTRUE(file.symlink(png_path, notebook_png_path))) {
    stop("Failed to link notebook figure: ", notebook_png_path, call. = FALSE)
  }
  notebook_png_path
}
