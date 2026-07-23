#!/usr/bin/env Rscript

# Targeted progenitor-high, proliferation-low proportion sensitivity analysis.
#
# Usage:
# ESPI_OVERWRITE=false \
# Rscript analysis/targeted-neurogenic-proportion/analysis.R
#
# Set ESPI_MODULE_SCORE_TABLE to override the upstream module-score table.

suppressPackageStartupMessages({
  library(here)
  here::i_am("analysis/targeted-neurogenic-proportion/analysis.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
  library(glmmTMB)
})

# ---- parameters ----

config <- publication_config()
bootstrap_replicates <- 1000L
progenitor_percentiles <- c(0.80, 0.90, 0.95)
proliferation_percentiles <- c(0.25, 0.50)
# ANALYSIS_OK[random-seeds]: fixed independent seeds make each gate bootstrap reproducible.
bootstrap_seeds <- c(6247L, 3518L, 9076L, 4661L, 7384L, 1859L)

primary_module_score_dir <- file.path(
  config$paths$degs,
  "mg_selected",
  "module_score_milo_da",
  "k_60__prop_0.04"
)
default_input_path <- file.path(primary_module_score_dir, "module_scores.tsv")
input_path <- Sys.getenv("ESPI_MODULE_SCORE_TABLE", unset = default_input_path)
output_dir <- here::here(
  "analysis",
  "targeted-neurogenic-proportion",
  "outputs"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_paths <- file.path(
  output_dir,
  c(
    "threshold_grid_results.tsv",
    "threshold_grid_sample_proportions.tsv",
    "bootstrap_null_statistics.tsv",
    "bootstrap_settings.tsv",
    "threshold_grid_sample_proportions.png",
    "threshold_grid_sample_proportions.pdf"
  )
)
assert_output_available(output_paths, config$overwrite)

if (!file.exists(input_path)) {
  cli::cli_abort(c(
    "The upstream module-score table does not exist.",
    "x" = "Missing path: {input_path}",
    "i" = paste(
      "Run scripts/module-score-milo-da.R first or set",
      "ESPI_MODULE_SCORE_TABLE."
    )
  ))
}

# ---- input ----

module_scores <- readr::read_tsv(input_path, show_col_types = FALSE)
required_columns <- c(
  "cell",
  "Mouse",
  "condition_label",
  "condition",
  "sample_id",
  "progenitor_score",
  "proliferation_score"
)
missing_columns <- base::setdiff(required_columns, colnames(module_scores))
if (length(missing_columns) > 0L) {
  cli::cli_abort(c(
    "The module-score table is missing required columns.",
    "x" = "Missing: {paste(missing_columns, collapse = ', ')}"
  ))
}
if (anyNA(module_scores[required_columns])) {
  cli::cli_abort("Required module-score columns contain missing values.")
}
if (anyDuplicated(module_scores$cell)) {
  cli::cli_abort("Cell identifiers must be unique.")
}
if (!base::setequal(unique(module_scores$condition), c("control", "estim"))) {
  cli::cli_abort("Condition must contain exactly control and estim.")
}

module_scores <- module_scores |>
  dplyr::mutate(
    Mouse = as.character(Mouse),
    condition = factor(condition, levels = c("control", "estim"))
  )

control_scores <- module_scores |> dplyr::filter(condition == "control")

threshold_grid <- tidyr::crossing(
  progenitor_percentile = progenitor_percentiles,
  proliferation_percentile = proliferation_percentiles
) |>
  dplyr::mutate(
    gate_index = dplyr::row_number(),
    bootstrap_seed = bootstrap_seeds,
    gate = paste0(
      "progenitor_p",
      100 * progenitor_percentile,
      "__proliferation_p",
      100 * proliferation_percentile
    )
  )

# ---- helpers ----

# ANALYSIS_OK[R026]: local helper derives equal-sample-weighted control quantiles.
weighted_empirical_quantile <- function(x, sample_id, probability) {
  sample_sizes <- table(sample_id)
  weights <- 1 / (length(sample_sizes) * as.numeric(sample_sizes[sample_id]))
  value_order <- order(x)
  sorted_values <- x[value_order]
  cumulative_weight <- cumsum(weights[value_order])
  sorted_values[which(cumulative_weight >= probability)[[1L]]]
}

# ANALYSIS_OK[R026]: local helper constructs one prespecified gate and its sample counts.
build_gate_data <- function(
  gate,
  gate_index,
  bootstrap_seed,
  progenitor_percentile,
  proliferation_percentile
) {
  progenitor_threshold <- weighted_empirical_quantile(
    control_scores$progenitor_score,
    control_scores$sample_id,
    progenitor_percentile
  )
  proliferation_threshold <- weighted_empirical_quantile(
    control_scores$proliferation_score,
    control_scores$sample_id,
    proliferation_percentile
  )

  sample_proportions <- module_scores |>
    dplyr::mutate(
      targeted = progenitor_score > progenitor_threshold &
        proliferation_score < proliferation_threshold
    ) |>
    dplyr::group_by(sample_id, Mouse, condition_label, condition) |>
    dplyr::summarise(
      responder_n = sum(targeted),
      sample_total = dplyr::n(),
      proportion = responder_n / sample_total,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      gate,
      gate_index,
      bootstrap_seed,
      progenitor_percentile,
      proliferation_percentile,
      progenitor_threshold,
      proliferation_threshold,
      .before = 1
    )

  model_data <- sample_proportions |>
    dplyr::mutate(
      condition = factor(condition, levels = c("control", "estim")),
      response = cbind(responder_n, sample_total - responder_n)
    )

  list(
    gate = gate,
    gate_index = gate_index,
    bootstrap_seed = bootstrap_seed,
    progenitor_percentile = progenitor_percentile,
    proliferation_percentile = proliferation_percentile,
    progenitor_threshold = progenitor_threshold,
    proliferation_threshold = proliferation_threshold,
    sample_proportions = sample_proportions,
    model_data = model_data
  )
}

# ANALYSIS_OK[R026]: local helper retries sparse beta-binomial fits with BFGS.
fit_bootstrap_response <- function(response, null_fit, full_fit, model_data) {
  default_result <- tryCatch(
    {
      null_refit <- suppressWarnings(glmmTMB::refit(null_fit, response))
      full_refit <- suppressWarnings(glmmTMB::refit(full_fit, response))
      if (isTRUE(null_refit$sdr$pdHess) && isTRUE(full_refit$sdr$pdHess)) {
        return(c(
          likelihood_ratio = max(
            0,
            2 *
              (as.numeric(stats::logLik(full_refit)) -
                as.numeric(stats::logLik(null_refit)))
          ),
          retried = 0,
          success = 1
        ))
      }
      NULL
    },
    error = function(error) NULL
  )
  if (!is.null(default_result)) {
    return(default_result)
  }

  tryCatch(
    {
      simulated_data <- model_data
      simulated_data$response <- response
      retry_control <- glmmTMB::glmmTMBControl(
        optimizer = stats::optim,
        optArgs = list(method = "BFGS")
      )
      # ANALYSIS_OK[contrast-definition]: one condition coefficient is the prespecified targeted contrast.
      null_refit <- suppressWarnings(glmmTMB::glmmTMB(
        response ~ 1,
        family = glmmTMB::betabinomial(link = "logit"),
        data = simulated_data,
        control = retry_control
      ))
      full_refit <- suppressWarnings(glmmTMB::glmmTMB(
        response ~ condition,
        family = glmmTMB::betabinomial(link = "logit"),
        data = simulated_data,
        control = retry_control
      ))
      if (isTRUE(null_refit$sdr$pdHess) && isTRUE(full_refit$sdr$pdHess)) {
        c(
          likelihood_ratio = max(
            0,
            2 *
              (as.numeric(stats::logLik(full_refit)) -
                as.numeric(stats::logLik(null_refit)))
          ),
          retried = 1,
          success = 1
        )
      } else {
        c(likelihood_ratio = NA_real_, retried = 1, success = 0)
      }
    },
    error = function(error) {
      c(likelihood_ratio = NA_real_, retried = 1, success = 0)
    }
  )
}

# ANALYSIS_OK[R026]: local helper runs and records one gate-specific null bootstrap.
bootstrap_gate <- function(gate_data) {
  model_data <- gate_data$model_data
  # ANALYSIS_OK[contrast-definition]: one condition coefficient is the prespecified targeted contrast.
  null_fit <- glmmTMB::glmmTMB(
    response ~ 1,
    family = glmmTMB::betabinomial(link = "logit"),
    data = model_data
  )
  full_fit <- glmmTMB::glmmTMB(
    response ~ condition,
    family = glmmTMB::betabinomial(link = "logit"),
    data = model_data
  )
  if (!isTRUE(null_fit$sdr$pdHess) || !isTRUE(full_fit$sdr$pdHess)) {
    cli::cli_abort("Observed beta-binomial fit has a non-positive Hessian.")
  }

  observed_lr <- max(
    0,
    2 *
      (as.numeric(stats::logLik(full_fit)) -
        as.numeric(stats::logLik(null_fit)))
  )
  coefficient <- summary(full_fit)$coefficients$cond["conditionestim", ]

  set.seed(gate_data$bootstrap_seed)
  simulations <- stats::simulate(
    null_fit,
    nsim = bootstrap_replicates,
    seed = gate_data$bootstrap_seed
  )
  simulation_results <- vapply(
    simulations,
    fit_bootstrap_response,
    numeric(3L),
    null_fit = null_fit,
    full_fit = full_fit,
    model_data = model_data
  )

  null_statistics <- tibble::tibble(
    gate = gate_data$gate,
    gate_index = gate_data$gate_index,
    bootstrap_seed = gate_data$bootstrap_seed,
    simulation = seq_len(bootstrap_replicates),
    likelihood_ratio = simulation_results["likelihood_ratio", ],
    retried = as.logical(simulation_results["retried", ]),
    success = as.logical(simulation_results["success", ])
  )
  valid_lr <- null_statistics$likelihood_ratio[
    null_statistics$success & !is.na(null_statistics$likelihood_ratio)
  ]
  exceedances <- sum(valid_lr >= observed_lr - 1e-12)
  failed_simulations <- bootstrap_replicates - length(valid_lr)
  bootstrap_p <- (exceedances + 1) / (length(valid_lr) + 1)

  condition_summary <- model_data |>
    dplyr::group_by(condition) |>
    dplyr::summarise(mean_proportion = mean(proportion), .groups = "drop")
  control_mean <- condition_summary$mean_proportion[
    condition_summary$condition == "control"
  ]
  estim_mean <- condition_summary$mean_proportion[
    condition_summary$condition == "estim"
  ]

  result <- tibble::tibble(
    gate = gate_data$gate,
    gate_index = gate_data$gate_index,
    bootstrap_seed = gate_data$bootstrap_seed,
    progenitor_percentile = gate_data$progenitor_percentile,
    proliferation_percentile = gate_data$proliferation_percentile,
    progenitor_threshold = gate_data$progenitor_threshold,
    proliferation_threshold = gate_data$proliferation_threshold,
    control_mean_percent = 100 * control_mean,
    estim_mean_percent = 100 * estim_mean,
    difference_percentage_points = 100 * (estim_mean - control_mean),
    ratio_of_mean_proportions = estim_mean / control_mean,
    odds_ratio = exp(coefficient[["Estimate"]]),
    odds_ratio_wald_low = exp(
      coefficient[["Estimate"]] - 1.96 * coefficient[["Std. Error"]]
    ),
    odds_ratio_wald_high = exp(
      coefficient[["Estimate"]] + 1.96 * coefficient[["Std. Error"]]
    ),
    observed_likelihood_ratio = observed_lr,
    requested_simulations = bootstrap_replicates,
    successful_simulations = length(valid_lr),
    retried_simulations = sum(null_statistics$retried),
    failed_simulations,
    exceedances,
    bootstrap_p,
    bootstrap_p_failure_low = (exceedances + 1) / (bootstrap_replicates + 1),
    bootstrap_p_failure_high = (exceedances + failed_simulations + 1) /
      (bootstrap_replicates + 1),
    monte_carlo_se = sqrt(bootstrap_p * (1 - bootstrap_p) / length(valid_lr))
  )

  list(result = result, null_statistics = null_statistics)
}

# ---- analysis ----

gate_data <- purrr::pmap(threshold_grid, build_gate_data)
bootstrap_outputs <- purrr::map(gate_data, bootstrap_gate)

threshold_results <- bootstrap_outputs |>
  purrr::map("result") |>
  purrr::list_rbind() |>
  dplyr::mutate(bootstrap_q_bh = stats::p.adjust(bootstrap_p, method = "BH"))

sample_proportions <- gate_data |>
  purrr::map("sample_proportions") |>
  purrr::list_rbind()

null_statistics <- bootstrap_outputs |>
  purrr::map("null_statistics") |>
  purrr::list_rbind()

paired_mice <- sample_proportions |>
  dplyr::distinct(Mouse, condition) |>
  dplyr::count(Mouse, name = "n_conditions") |>
  dplyr::filter(n_conditions == 2L) |>
  dplyr::pull(Mouse)

plot_data <- sample_proportions |>
  dplyr::mutate(
    condition = factor(condition, levels = c("control", "estim")),
    progenitor_gate = factor(
      progenitor_percentile,
      levels = progenitor_percentiles,
      labels = paste0("Progenitor > p", 100 * progenitor_percentiles)
    ),
    proliferation_gate = factor(
      proliferation_percentile,
      levels = proliferation_percentiles,
      labels = paste0("Proliferation < p", 100 * proliferation_percentiles)
    )
  )

plot_sample_proportions <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = condition, y = proportion)
) +
  ggplot2::geom_line(
    data = \(data) dplyr::filter(data, Mouse %in% paired_mice),
    ggplot2::aes(group = Mouse),
    linewidth = 0.6
  ) +
  ggplot2::geom_point(size = 2.5) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = Mouse),
    direction = "y",
    seed = config$seed,
    min.segment.length = 0,
    box.padding = 0.15,
    point.padding = 0.1,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(proliferation_gate),
    cols = ggplot2::vars(progenitor_gate)
  ) +
  ggplot2::scale_x_discrete(
    labels = c(control = "p27CKO", estim = "p27CKO + EStim")
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.04, 0.16))
  ) +
  ggplot2::labs(x = NULL, y = "Cells meeting both score thresholds") +
  theme_stone() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
    panel.spacing = grid::unit(1, "lines")
  )

settings <- tibble::tibble(
  parameter = c(
    "default_input_path",
    "input_path_used",
    "input_override_used",
    "input_sha256",
    "output_dir",
    "analysis_unit",
    "model",
    "link",
    "control_quantile_weighting",
    "progenitor_percentiles",
    "proliferation_percentiles",
    "bootstrap_replicates_per_gate",
    "bootstrap_seeds",
    "failed_fit_retry_optimizer",
    "bootstrap_p_correction",
    "glmmTMB_version",
    "R_version"
  ),
  value = c(
    default_input_path,
    input_path,
    as.character(!identical(input_path, default_input_path)),
    digest::digest(file = input_path, algo = "sha256"),
    file.path("analysis", "targeted-neurogenic-proportion", "outputs"),
    "Mouse x Condition sample",
    "beta-binomial condition-only likelihood-ratio test",
    "logit",
    "equal weight per control sample",
    paste(progenitor_percentiles, collapse = ","),
    paste(proliferation_percentiles, collapse = ","),
    as.character(bootstrap_replicates),
    paste(bootstrap_seeds, collapse = ","),
    "BFGS",
    "BH across six prespecified gates",
    as.character(utils::packageVersion("glmmTMB")),
    R.version.string
  )
)

# ---- outputs ----

readr::write_tsv(
  threshold_results,
  file.path(output_dir, "threshold_grid_results.tsv")
)
readr::write_tsv(
  sample_proportions,
  file.path(output_dir, "threshold_grid_sample_proportions.tsv")
)
readr::write_tsv(
  null_statistics,
  file.path(output_dir, "bootstrap_null_statistics.tsv")
)
readr::write_tsv(settings, file.path(output_dir, "bootstrap_settings.tsv"))
ggplot2::ggsave(
  file.path(output_dir, "threshold_grid_sample_proportions.png"),
  plot_sample_proportions,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 300
)
ggplot2::ggsave(
  file.path(output_dir, "threshold_grid_sample_proportions.pdf"),
  plot_sample_proportions,
  width = 10,
  height = 7.5,
  units = "in"
)
