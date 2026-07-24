# Scoring cells along reactive/reprogramming axis and proliferation axis

devtools::load_all()

library(here)
library(tidyverse)
library(Seurat)
library(ggview)

here::i_am("scripts/reprog_scoring.R")

# --- helpers
robust_z_score <- function(x) {
  (x - median(x)) / mad(x)
}

# ANALYSIS_OK[script-helper]: called by observed and bootstrap decompositions below.
decompose_module_scores <- function(sample_cluster_scores, sample_weights) {
  weighted_scores <- sample_weights |>
    dplyr::left_join(
      sample_cluster_scores,
      by = c("Sample", "condition_role"),
      relationship = "many-to-many"
    )
  combinations_per_sample <- dplyr::n_distinct(sample_cluster_scores$cluster) *
    dplyr::n_distinct(sample_cluster_scores$module)
  stopifnot(
    nrow(weighted_scores) == nrow(sample_weights) * combinations_per_sample
  )

  weighted_scores |>
    dplyr::group_by(draw, condition_role, cluster, module) |>
    dplyr::summarize(
      proportion = sum(sample_weight * cluster_proportion),
      score_contribution = sum(sample_weight * cluster_score_contribution),
      .groups = "drop"
    ) |>
    dplyr::mutate(within_mean = score_contribution / proportion) |>
    dplyr::select(-score_contribution) |>
    tidyr::pivot_wider(
      names_from = condition_role,
      values_from = c(proportion, within_mean),
      names_glue = "{.value}_{condition_role}"
    ) |>
    dplyr::mutate(
      composition = 0.5 *
        (proportion_estim - proportion_control) *
        (within_mean_estim + within_mean_control),
      within_state = 0.5 *
        (proportion_estim + proportion_control) *
        (within_mean_estim - within_mean_control),
      total = composition + within_state
    )
}

# ---- parameters ----

config <- publication_config()
input_path <- config$selected$mg$path
condition_col <- config$conditions$column
cluster_col <- config$selected$mg$column
control_label <- config$conditions$control
estim_label <- config$conditions$estim
seed <- config$seed
decomposition_draws <- 5000L
decomposition_tolerance <- 1e-10

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
    cluster = .data[[cluster_col]],
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
mod_score_summarized <- mod_scores |>
  group_by(Sample, Mouse, Condition) |>
  summarize(
    across(
      c(neuro_progen_score, prolif_score, mg_score, cone_score),
      mean,
      .names = "{.col}"
    ),
    .groups = "drop"
  )

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

module_labels <- c(
  neuro_progen_score = "Neural progenitor",
  prolif_score = "Proliferation",
  mg_score = "Müller glia",
  cone_score = "Cone bipolar"
)
decomposition_score_columns <- names(module_labels)

unexpected_conditions <- setdiff(
  unique(mod_scores$Condition),
  c(control_label, estim_label)
)
if (length(unexpected_conditions) > 0L) {
  stop("Unexpected condition labels in module-score data.", call. = FALSE)
}

sample_cell_counts <- mod_scores |>
  dplyr::count(Sample, Condition, name = "sample_cells")

decomposition_sample_cluster_scores <- mod_scores |>
  dplyr::select(
    Sample,
    Condition,
    cluster,
    dplyr::all_of(decomposition_score_columns)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(decomposition_score_columns),
    names_to = "module",
    values_to = "score"
  ) |>
  dplyr::group_by(Sample, Condition, cluster, module) |>
  dplyr::summarize(
    cluster_cells = dplyr::n(),
    score_sum = sum(score),
    .groups = "drop"
  ) |>
  tidyr::complete(
    tidyr::nesting(Sample, Condition),
    cluster,
    module,
    fill = list(cluster_cells = 0L, score_sum = 0)
  )

sample_cluster_rows <- nrow(decomposition_sample_cluster_scores)
decomposition_sample_cluster_scores <- decomposition_sample_cluster_scores |>
  dplyr::left_join(sample_cell_counts, by = c("Sample", "Condition"))
stopifnot(nrow(decomposition_sample_cluster_scores) == sample_cluster_rows)

decomposition_sample_cluster_scores <- decomposition_sample_cluster_scores |>
  dplyr::mutate(
    condition_role = dplyr::case_when(
      Condition == control_label ~ "control",
      Condition == estim_label ~ "estim"
    ),
    cluster_proportion = cluster_cells / sample_cells,
    cluster_score_contribution = score_sum / sample_cells
  ) |>
  dplyr::select(
    Sample,
    condition_role,
    cluster,
    module,
    cluster_proportion,
    cluster_score_contribution
  )

decomposition_sample_table <- decomposition_sample_cluster_scores |>
  dplyr::distinct(Sample, condition_role)

decomposition_observed_weights <- decomposition_sample_table |>
  dplyr::group_by(condition_role) |>
  dplyr::mutate(draw = 0L, sample_weight = 1 / dplyr::n()) |>
  dplyr::ungroup()

decomposition_by_cluster <- decompose_module_scores(
  decomposition_sample_cluster_scores,
  decomposition_observed_weights
)

decomposition_by_module <- decomposition_by_cluster |>
  dplyr::group_by(module) |>
  dplyr::summarize(
    composition = sum(composition),
    within_state = sum(within_state),
    total = sum(total),
    .groups = "drop"
  )

decomposition_direct_contrasts <- mod_scores |>
  dplyr::select(
    Sample,
    Condition,
    dplyr::all_of(decomposition_score_columns)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(decomposition_score_columns),
    names_to = "module",
    values_to = "score"
  ) |>
  dplyr::group_by(Sample, Condition, module) |>
  dplyr::summarize(score = mean(score), .groups = "drop") |>
  dplyr::group_by(Condition, module) |>
  dplyr::summarize(score = mean(score), .groups = "drop") |>
  dplyr::mutate(
    condition_role = dplyr::if_else(
      Condition == estim_label,
      "estim",
      "control"
    )
  ) |>
  dplyr::select(-Condition) |>
  tidyr::pivot_wider(names_from = condition_role, values_from = score) |>
  dplyr::mutate(direct_total = estim - control)

module_decomposition_rows <- nrow(decomposition_by_module)
decomposition_by_module <- decomposition_by_module |>
  dplyr::left_join(
    decomposition_direct_contrasts |> dplyr::select(module, direct_total),
    by = "module"
  )
stopifnot(nrow(decomposition_by_module) == module_decomposition_rows)

decomposition_by_module <- decomposition_by_module |>
  dplyr::mutate(reconstruction_error = total - direct_total)

if (
  max(abs(decomposition_by_module$reconstruction_error)) >
    decomposition_tolerance
) {
  stop(
    "Composition-state decomposition failed to reconstruct totals.",
    call. = FALSE
  )
}

set.seed(seed)
decomposition_bootstrap_weights <- decomposition_sample_table |>
  dplyr::group_by(condition_role) |>
  dplyr::group_modify(\(data, key) {
    # ANALYSIS_OK[random-seed-only]: RNG generates Bayesian-bootstrap weights.
    raw_weights <- matrix(
      stats::rexp(decomposition_draws * nrow(data)),
      nrow = decomposition_draws,
      ncol = nrow(data)
    )
    weights <- raw_weights / rowSums(raw_weights)

    tibble::tibble(
      draw = rep(seq_len(decomposition_draws), each = nrow(data)),
      Sample = rep(data$Sample, times = decomposition_draws),
      sample_weight = as.vector(t(weights))
    )
  }) |>
  dplyr::ungroup()

decomposition_bootstrap_by_cluster <- decompose_module_scores(
  decomposition_sample_cluster_scores,
  decomposition_bootstrap_weights
)

decomposition_bootstrap_by_module <- decomposition_bootstrap_by_cluster |>
  dplyr::group_by(draw, module) |>
  dplyr::summarize(
    composition = sum(composition),
    within_state = sum(within_state),
    total = sum(total),
    .groups = "drop"
  )

decomposition_summary <- decomposition_bootstrap_by_module |>
  tidyr::pivot_longer(
    cols = c(composition, within_state, total),
    names_to = "component",
    values_to = "value"
  ) |>
  dplyr::group_by(module, component) |>
  dplyr::summarize(
    lower_95 = stats::quantile(value, 0.025),
    upper_95 = stats::quantile(value, 0.975),
    probability_positive = mean(value > 0),
    .groups = "drop"
  )

decomposition_estimates <- decomposition_by_module |>
  dplyr::select(module, composition, within_state, total) |>
  tidyr::pivot_longer(
    cols = c(composition, within_state, total),
    names_to = "component",
    values_to = "estimate"
  )

decomposition_summary_rows <- nrow(decomposition_summary)
decomposition_summary <- decomposition_summary |>
  dplyr::left_join(decomposition_estimates, by = c("module", "component"))
stopifnot(nrow(decomposition_summary) == decomposition_summary_rows)

decomposition_summary <- decomposition_summary |>
  dplyr::mutate(
    module = factor(
      module,
      levels = names(module_labels),
      labels = unname(module_labels)
    ),
    component = factor(
      component,
      levels = c("composition", "within_state", "total"),
      labels = c("Composition", "Within-state", "Total")
    )
  )

decomposition_summary

decomposition_plot <- decomposition_summary |>
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = module)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0
  ) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(ggplot2::vars(component), ncol = 1) +
  ggplot2::labs(
    x = "E-Stim minus control score difference",
    y = NULL,
    title = "Composition-versus-state decomposition",
    subtitle = "Intervals use a Bayesian bootstrap over samples"
  ) +
  theme_stone()

decomposition_plot

decomposition_cluster_plot <- decomposition_by_cluster |>
  dplyr::select(module, cluster, composition, within_state) |>
  tidyr::pivot_longer(
    cols = c(composition, within_state),
    names_to = "component",
    values_to = "contribution"
  ) |>
  dplyr::mutate(
    module = factor(
      module,
      levels = names(module_labels),
      labels = unname(module_labels)
    ),
    component = factor(
      component,
      levels = c("composition", "within_state"),
      labels = c("Composition", "Within-state")
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = cluster, y = module, fill = contribution)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(ggplot2::vars(component)) +
  ggplot2::scale_fill_gradient2(midpoint = 0) +
  ggplot2::labs(
    x = "Cluster",
    y = NULL,
    fill = "Contribution",
    title = "Cluster contributions to module-score differences"
  ) +
  theme_stone()

decomposition_cluster_plot

axis_pca <- prcomp(
  mod_scores[c("mg_score", "neuro_progen_score", "prolif_score", "cone_score")],
  center = TRUE,
  scale. = TRUE
)

autoplot(axis_pca, data = mod_scores, color = "Condition")

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


fit_scores <- lm(
  cbind(mg_score, neuro_progen_score, prolif_score, cone_score) ~ Condition,
  data = mod_scores
)

fit_scores_persample <- lm(
  cbind(mg_score, neuro_progen_score, prolif_score, cone_score) ~ Condition,
  data = mod_score_summarized
)

score_model_prior <- c(
  brms::set_prior("normal(0, 5)", class = "b", resp = "mgscore"),
  brms::set_prior("normal(0, 5)", class = "b", resp = "neuroprogenscore"),
  brms::set_prior("normal(0, 5)", class = "b", resp = "prolifscore"),
  brms::set_prior("normal(0, 5)", class = "b", resp = "conescore")
)

fit_scores_randomsample <- brms::brm(
  brms::bf(
    brms::mvbind(
      mg_score,
      neuro_progen_score,
      prolif_score,
      cone_score
    ) ~ Condition + (1 | Sample)
  ) +
    brms::set_rescor(TRUE),
  data = mod_scores,
  prior = score_model_prior,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed
)

condition_effects <- brms::fixef(
  fit_scores_randomsample,
  probs = c(0.025, 0.975)
) |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  dplyr::filter(stringr::str_detect(term, "_Condition")) |>
  dplyr::transmute(
    score = stringr::str_remove(term, "_Condition.*$"),
    estimate = Estimate,
    posterior_sd = Est.Error,
    lower_95 = Q2.5,
    upper_95 = Q97.5,
    probability_positive = purrr::map_dbl(term, \(x) {
      mean(
        posterior::as_draws_df(fit_scores_randomsample)[[paste0("b_", x)]] > 0
      )
    })
  ) |>
  dplyr::mutate(
    score = dplyr::recode(
      score,
      mgscore = "MG",
      neuroprogenscore = "Neural progenitor",
      prolifscore = "Proliferation",
      conescore = "Cone"
    )
  )
condition_effects

condition_effects |>
  dplyr::mutate(score = forcats::fct_reorder(score, estimate)) |>
  ggplot2::ggplot(ggplot2::aes(y = score, x = estimate)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0
  ) +
  ggplot2::geom_point() +
  ggplot2::labs(
    x = "E-Stim effect (standard deviations)",
    y = NULL,
    title = "Posterior condition effects",
    subtitle = "Posterior means and 95% credible intervals"
  ) +
  theme_stone()
