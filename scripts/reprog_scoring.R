# Scoring cells along reactive/reprogramming axis and proliferation axis

devtools::load_all()

library(here)
library(tidyverse)
library(Seurat)

here::i_am("scripts/reprog_scoring.R")

# --- helpers
robust_z_score <- function(x) {
  (x - median(x)) / mad(x)
}

# ---- parameters ----

config <- publication_config()
input_path <- config$selected$mg$path
condition_col <- config$conditions$column
control_label <- config$conditions$control
estim_label <- config$conditions$estim
seed <- config$seed

# Load sobj
sobj <- readRDS(input_path)

marker_table <- stack(cell_type_marker_genes) |>
  tibble::as_tibble() |>
  dplyr::rename(gene = values, cell_type = ind) |>
  dplyr::mutate(
    cell_type = as.character(cell_type),
    cell_type_label = base::unname(cell_type_marker_labels[cell_type])
  ) |>
  dplyr::bind_rows(tibble::tibble(
    gene = "Cdkn1b",
    cell_type = "cdkn1b_standalone",
    cell_type_label = "Cdkn1b"
  ))

sobj <- AddModuleScore(
  object = sobj,
  features = list(
    neuro_progen = marker_table |>
      dplyr::filter(cell_type == "neurogenic_progenitor") |>
      dplyr::pull(gene),
    prolif = marker_table |>
      dplyr::filter(cell_type == "proliferative") |>
      dplyr::pull(gene),
    mg = marker_table |>
      dplyr::filter(cell_type == "muller_glia") |>
      dplyr::pull(gene),
    activated_mg = marker_table |>
      dplyr::filter(cell_type == "activated_muller_glia") |>
      dplyr::pull(gene),
    cone = marker_table |>
      dplyr::filter(cell_type == "cone_bipolar") |>
      dplyr::pull(gene)
  ),
  name = c("neuro_progen", "prolif", "mg", "activated_mg", "cone"),
  assay = "RNA",
  search = FALSE,
  seed = seed
)

mod_scores <- sobj[[]] |>
  tibble::rownames_to_column("cell") |>
  dplyr::mutate(
    cell,
    Sample,
    Mouse,
    Condition,
    neuro_progen_score = neuro_progen1,
    prolif_score = prolif2,
    mg_score = mg3,
    activated_mg_score = activated_mg4,
    cone_score = cone5,
    .keep = "none"
  )

mod_scores <- mod_scores |>
  dplyr::mutate(mg_score_raw = mg_score, cone_score_raw = cone_score) |>
  dplyr::mutate(dplyr::across(ends_with("score"), \(x) as.numeric(scale(x))))

plot_mod_scores <- mod_scores |>
  pivot_longer(
    cols = c(
      neuro_progen_score,
      prolif_score,
      mg_score,
      activated_mg_score,
      cone_score
    ),
    names_to = "module",
    values_to = "score"
  ) |>
  ggplot(aes(x = Condition, y = score, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.2) +
  facet_wrap(~module, scales = "free_y") +
  theme_bw()

plot_mod_scores

scorelist <- grep("score$", colnames(mod_scores), value = TRUE)

neuron_prolif_scores <- mod_scores |>
  dplyr::mutate(neuro_mg = neuro_progen_score + cone_score - mg_score)


plot_mg_v_activated <- mod_scores |>
  ggplot(aes(
    x = mg_score,
    y = activated_mg_score,
    color = Condition,
    fill = Condition
  )) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth() +
  scale_color_manual(
    values = c(
      "p27CKO" = "gray35",
      "p27CKO +EStim" = config$palettes$dotplot[[2]]
    ),
    aesthetics = c("color", "fill")
  ) +
  theme_bw()

plot_mg_v_activated

mg_score_limits <- range(mod_scores$mg_score_raw)
cone_score_limits <- range(mod_scores$cone_score_raw)

plot_mg_cone_density <- mod_scores |>
  dplyr::mutate(
    Condition = factor(Condition, levels = c(control_label, estim_label))
  ) |>
  ggplot(aes(x = mg_score_raw, y = cone_score_raw)) +
  geom_hex() +
  facet_wrap(~Condition) +
  coord_cartesian(
    xlim = mg_score_limits,
    ylim = cone_score_limits,
    expand = FALSE
  ) +
  scale_fill_viridis_c() +
  labs(x = "Müller glia module score", y = "Cone bipolar module score") +
  theme_stone()

plot_mg_cone_density

plot_progen_cone <- mod_scores |>
  ggplot(aes(x = neuro_progen_score, y = cone_score)) +
  facet_wrap(~Condition) +
  geom_hex() +
  scale_fill_viridis_c() +
  theme_bw()

plot_progen_cone

plot_mg_v_prolif <- mod_scores |>
  ggplot(aes(
    x = mg_score,
    y = prolif_score,
    color = Condition,
    fill = Condition
  )) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth() +
  scale_color_manual(
    values = c(
      "p27CKO" = "gray35",
      "p27CKO +EStim" = config$palettes$dotplot[[2]]
    ),
    aesthetics = c("color", "fill")
  ) +
  theme_bw()

plot_mg_v_prolif

plot_mg_progen <- mod_scores |>
  ggplot(aes(
    x = mg_score,
    y = neuro_progen_score,
    color = Condition,
    fill = Condition
  )) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth() +
  scale_color_manual(
    values = c(
      "p27CKO" = "gray35",
      "p27CKO +EStim" = config$palettes$dotplot[[2]]
    ),
    aesthetics = c("color", "fill")
  ) +
  theme_bw()

plot_mg_progen

# Now per sample for stats
mods_samples_stats <- mod_scores |>
  mutate(
    Condition = factor(Condition, levels = c(control_label, estim_label))
  ) |>
  group_by(Sample, Mouse, Condition) |>
  summarize(
    across(
      c(neuro_progen_score, prolif_score, mg_score, activated_mg_score),
      mean,
      .names = "{.col}"
    ),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(neuro_progen_score, prolif_score, mg_score, activated_mg_score),
    names_to = "module",
    values_to = "score"
  ) |>
  nest_by(module) |>
  mutate(
    t_test = list(t.test(score ~ Condition, data = data)),
    tidy_t_test = list(broom::tidy(t_test))
  ) |>
  unnest(tidy_t_test) |>
  dplyr::select(
    module,
    estimate1,
    estimate2,
    statistic,
    p.value,
    conf.low,
    conf.high
  ) |>
  dplyr::mutate(
    module = factor(
      module,
      levels = c(
        "neuro_progen_score",
        "prolif_score",
        "mg_score",
        "activated_mg_score"
      )
    )
  )

plot_mods_samples <- mod_scores |>
  pivot_longer(
    cols = c(neuro_progen_score, prolif_score, mg_score, activated_mg_score),
    names_to = "module",
    values_to = "score"
  ) |>
  group_by(Sample, Mouse, Condition, module) |>
  summarize(score = mean(score), .groups = "drop") |>
  ggplot(aes(x = Condition, y = score, color = Condition)) +
  facet_wrap(~module, scales = "free_y") +
  geom_jitter(alpha = 0.5, width = 0.2, height = 0) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.25, ) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  scale_color_manual(
    values = c(
      "p27CKO" = "gray35",
      "p27CKO +EStim" = config$palettes$dotplot[[2]]
    ),
    aesthetics = c("color", "fill")
  ) +
  theme_stone() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

plot_mods_samples + canvas(6, 4)


axis_pca <- prcomp(
  mod_scores[c("mg_score", "neuro_progen_score", "prolif_score", "cone_score")],
  center = TRUE,
  scale. = TRUE
)

axis_pca_points <- axis_pca$x |> as.data.frame() |> bind_cols(mod_scores)

mod_pca_plot <- axis_pca_points |>
  ggplot(aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(
    values = c(
      "p27CKO" = "gray35",
      "p27CKO +EStim" = config$palettes$dotplot[[2]]
    ),
    aesthetics = c("color", "fill")
  ) +
  theme_bw()

mod_pca_plot

mod_grid <- mod_scores[c(
  "mg_score",
  "activated_mg_score",
  "neuro_progen_score",
  "prolif_score"
)] |>
  summarize(across(
    everything(),
    .fns = list(min = min, median = median, max = max),
    .names = "{.col}.{.fn}"
  )) |>
  pivot_longer(
    cols = everything(),
    names_to = c("module", "stat"),
    names_sep = "\\.",
    values_to = "value"
  ) |>
  dplyr::filter(stat != "median") |>
  pivot_wider(names_from = module, values_from = value)
