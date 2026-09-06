test_that("sequential marker heatmaps preserve block and cell order", {
  withr::local_seed(6381)
  genes <- unique(unlist(cell_type_marker_genes, use.names = FALSE))
  counts <- matrix(
    sample.int(20L, length(genes) * 12L, replace = TRUE),
    nrow = length(genes),
    dimnames = list(genes, paste0("cell", seq_len(12L)))
  )
  sobj <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)
  )
  sobj$cluster <- rep(c("10", "2", "0", "1"), 3L)
  output_dir <- withr::local_tempdir()

  original_draw <- ComplexHeatmap::draw
  rendered <- NULL
  local_mocked_bindings(
    draw = function(object, ...) {
      drawn <- original_draw(object, ...)
      if (inherits(object, "Heatmap")) {
        rendered <<- drawn
      }
      invisible(drawn)
    },
    .package = "ComplexHeatmap"
  )

  modes <- list(
    sequential = list(sequential_clusters = TRUE),
    sequential_override = list(
      sequential_clusters = TRUE,
      cluster_cells = TRUE
    ),
    cluster_means = list(),
    cluster_cells = list(cluster_cells = TRUE)
  )
  for (mode in names(modes)) {
    paths <- do.call(
      write_curated_marker_heatmap,
      c(
        list(
          sobj = sobj,
          cluster_column = "cluster",
          layer = "counts",
          output_stem = file.path(output_dir, mode),
          width = 10,
          height = 9
        ),
        modes[[mode]]
      )
    )
    expect_equal(unname(file.exists(paths)), c(TRUE, TRUE))
    heatmap <- rendered@ht_list[["Row z-score"]]
    expect_length(heatmap@column_title, 0L)
    expect_equal(
      rownames(heatmap@matrix),
      unname(unlist(cell_type_marker_genes))
    )
    expect_identical(heatmap@column_dend_param$cluster, mode == "cluster_cells")
    expect_identical(
      heatmap@column_dend_param$cluster_slices,
      mode == "cluster_cells"
    )

    if (startsWith(mode, "sequential")) {
      order <- ComplexHeatmap::column_order(rendered)
      expect_identical(names(order), paste("Cluster", c("0", "1", "2", "10")))
      expect_equal(
        unname(order),
        lapply(c("0", "1", "2", "10"), function(id) {
          unname(which(sobj$cluster == id))
        })
      )
      expect_identical(names(heatmap@top_annotation), "Cluster")
    } else if (mode == "cluster_means") {
      expect_identical(
        names(heatmap@top_annotation),
        c("Cluster_dendrogram", "Cluster")
      )
    }
  }
})
