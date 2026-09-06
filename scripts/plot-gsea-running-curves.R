#!/usr/bin/env Rscript

# Supplemental GSEA curves reconstructed from saved results; no enrichment rerun.
# Usage: source("scripts/plot-gsea-running-curves.R")

suppressPackageStartupMessages({
  here::i_am("scripts/plot-gsea-running-curves.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# ---- parameters ----

config <- publication_config()
overwrite <- FALSE # Change only with permission to replace these outputs.
width <- 6.5
height <- 4.5
dpi <- 300
term_ids <- c("GO:0034341", "GO:0006805", "GO:1902969")
term_slugs <- c(
  "interferon_response",
  "xenobiotic_metabolism",
  "mitotic_dna_replication"
)
enrichment_dir <- file.path(config$paths$enrichment, "mg_selected")
figure_dir <- file.path(config$paths$figures, "mg_selected")

# ---- saved rankings and gene sets ----

results <- readr::read_tsv(
  file.path(enrichment_dir, "go_bp_gsea.tsv"),
  show_col_types = FALSE
)
terms <- results[match(term_ids, results$ID), ]
stopifnot(
  !anyDuplicated(results$ID[results$ID %in% term_ids]),
  identical(terms$ID, term_ids),
  all(is.finite(terms$enrichmentScore)),
  all(is.finite(terms$NES)),
  all(is.finite(terms$p.adjust))
)
ranks <- readr::read_tsv(
  file.path(enrichment_dir, "go_bp_gsea_symbol_entrez_mapping.tsv"),
  col_types = readr::cols(ENTREZID = readr::col_character())
) |>
  dplyr::filter(selected_for_gsea) |>
  dplyr::arrange(dplyr::desc(stat)) |>
  dplyr::mutate(position = dplyr::row_number())
stopifnot(
  nrow(ranks) > 0L,
  !anyNA(ranks$ENTREZID),
  !anyDuplicated(ranks$ENTREZID),
  all(is.finite(ranks$stat))
)

# GOALL includes ancestor annotations, matching gseGO's GO membership.
membership <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = term_ids,
  keytype = "GOALL",
  columns = "ENTREZID"
) |>
  dplyr::distinct(GOALL, ENTREZID)

# ---- reconstruct and validate ----

curves <- purrr::map(term_ids, function(term_id) {
  saved <- dplyr::filter(terms, ID == term_id)
  hits <- ranks$ENTREZID %in% membership$ENTREZID[membership$GOALL == term_id]
  weights <- abs(ranks$stat) * hits # Original gseGO exponent = 1.
  stopifnot(
    sum(hits) == saved$setSize,
    sum(hits) < length(hits),
    sum(weights) > 0
  )
  score <- cumsum(weights / sum(weights) - (!hits) / sum(!hits))
  peak <- which.max(abs(score))
  leading_positions <- if (saved$enrichmentScore > 0) {
    ranks$position <= peak
  } else {
    ranks$position > peak
  }
  leading_ids <- ranks$ENTREZID[hits & leading_positions]

  # Allow numerical rounding, but reject changed rankings or GO annotations.
  stopifnot(
    abs(score[peak] - saved$enrichmentScore) < 1e-6,
    base::setequal(
      leading_ids,
      strsplit(saved$core_enrichment, "/", fixed = TRUE)[[1]]
    ),
    abs(tail(score, 1)) < 1e-8
  )
  tibble::tibble(
    ID = term_id,
    position = c(0L, ranks$position),
    hit = c(FALSE, hits),
    running_score = c(0, score)
  )
}) |>
  purrr::list_rbind()

# ---- plots and Box exports ----

score_limits <- range(curves$running_score)
output_stems <- file.path(
  figure_dir,
  paste0("mg_selected_gsea_running_", term_slugs)
)
# Check every destination before writing any of the six files.
output_files <- output_path(
  c(paste0(output_stems, ".png"), paste0(output_stems, ".pdf")),
  overwrite = overwrite
)

for (i in seq_along(term_ids)) {
  term <- terms[i, ]
  curve <- dplyr::filter(curves, ID == term_ids[[i]])
  plot <- ggplot(curve, aes(position, running_score)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_line(linewidth = 0.7) +
    geom_rug(data = dplyr::filter(curve, hit), sides = "b", linewidth = 0.3) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(
      limits = score_limits,
      expand = expansion(mult = c(0.06, 0.08))
    ) +
    labs(
      title = term$Description,
      subtitle = sprintf(
        "NES = %.2f | BH-adjusted p = %.2g | %d genes",
        term$NES,
        term$p.adjust,
        term$setSize
      ),
      x = "Gene rank by DESeq2 Wald statistic\nE-Stim-associated (left) | p27CKO-associated (right)",
      y = "Running enrichment score",
      caption = paste(
        term$ID,
        "| Ticks: gene-set members | MG-selected pseudobulk"
      )
    ) +
    theme_stone(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8, hjust = 0),
      axis.title.x = element_text(size = 10, margin = margin(t = 8))
    )

  # Box-only supplemental exports; do not copy into the manuscript notebook.
  for (extension in c(".png", ".pdf")) {
    ggplot2::ggsave(
      filename = output_path(
        paste0(output_stems[[i]], extension),
        overwrite = overwrite
      ),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  }
}
