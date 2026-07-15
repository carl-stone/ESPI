set positional-arguments

# List available recipes
_default:
    just --list

# Run the canonical pipeline from counts-qc, legacy, or an explicit RDS
[group: "Canonical pipeline"]
run source="counts-qc" overwrite="false": (_run-pipeline source overwrite "false")

# Print the canonical pipeline plan without changing files
[group: "Canonical pipeline"]
run-dry-run source="counts-qc" overwrite="false": (_run-pipeline source overwrite "true")

# Build arguments for the canonical pipeline runner
_run-pipeline source="counts-qc" overwrite="false" dry_run="false":
    #!/usr/bin/env bash
    set -euo pipefail
    case "{{overwrite}}" in
        false|true)
            ;;
        *)
            echo "overwrite must be exactly true or false" >&2
            exit 2
            ;;
    esac
    args=()
    case "{{source}}" in
        counts-qc|legacy)
            args+=(--input-source "{{source}}")
            ;;
        *)
            args+=(--input "{{source}}")
            ;;
    esac
    if [ "{{dry_run}}" = "true" ]; then
        args+=(--dry-run)
    fi
    if [ "{{overwrite}}" = "true" ]; then
        args+=(--overwrite)
    fi
    Rscript scripts/run-pipeline.R "${args[@]}"


# Load package code
[group: "Expert and maintenance"]
load:
    Rscript -e 'devtools::load_all(".", quiet = TRUE)'

# Update package documentation
[group: "Expert and maintenance"]
document:
    Rscript -e 'devtools::document()'

# Rebuild README.md from README.Rmd
[group: "Expert and maintenance"]
readme:
    Rscript -e 'devtools::build_readme()'

# Format R code with Air
[group: "Expert and maintenance"]
format *paths:
    #!/usr/bin/env bash
    set -euo pipefail
    if [ "$#" -eq 0 ]; then
        air format R scripts
    else
        air format {{paths}}
    fi

# Run scilintr over first-party analysis code
# scilintr 0.1.1 accepts one root per invocation, so keep this loop scoped.
[group: "Expert and maintenance"]
lint:
    #!/usr/bin/env bash
    set -euo pipefail
    for scope in R scripts data-raw tools notebook/sc_analysis.qmd config.local.example.R; do
        Rscript -e 'scope <- commandArgs(trailingOnly = TRUE)[1]; quit(status = scilintr::main(scope))' "$scope"
    done

# Run all preprocessing branches from legacy or counts-qc input
[group: "Expert and maintenance"]
preprocess input_source="legacy" input="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--input-source "{{input_source}}")
    if [ -n "{{input}}" ]; then
        args=(--input "{{input}}")
    fi
    Rscript scripts/03-preprocess-all.R "${args[@]}"

# Run one preprocessing branch from legacy, counts-qc, or an explicit input
[group: "Expert and maintenance"]
preprocess-one normalization="pflog" filter_cc="false" input_source="legacy" input="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--normalization "{{normalization}}" --input-source "{{input_source}}")
    if [ -n "{{input}}" ]; then
        args=(--normalization "{{normalization}}" --input "{{input}}")
    fi
    if [ "{{filter_cc}}" = "true" ]; then
        args+=(--filter-cell-cycle)
    fi
    Rscript scripts/03-preprocess.R "${args[@]}"

# Cluster every preprocessed object
[group: "Expert and maintenance"]
cluster elbow_n="20" input_dir="" extra_dims="" resolutions="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--elbow-n "{{elbow_n}}")
    if [ -n "{{input_dir}}" ]; then
        args+=(--input-dir "{{input_dir}}")
    fi
    if [ -n "{{extra_dims}}" ]; then
        args+=(--extra-dims "{{extra_dims}}")
    fi
    if [ -n "{{resolutions}}" ]; then
        args+=(--resolutions "{{resolutions}}")
    fi
    Rscript scripts/04-cluster-all.R "${args[@]}"

# Print cluster commands without running them
[group: "Expert and maintenance"]
cluster-dry-run elbow_n="20" input_dir="" extra_dims="" resolutions="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--dry-run --elbow-n "{{elbow_n}}")
    if [ -n "{{input_dir}}" ]; then
        args+=(--input-dir "{{input_dir}}")
    fi
    if [ -n "{{extra_dims}}" ]; then
        args+=(--extra-dims "{{extra_dims}}")
    fi
    if [ -n "{{resolutions}}" ]; then
        args+=(--resolutions "{{resolutions}}")
    fi
    Rscript scripts/04-cluster-all.R "${args[@]}"

# Cluster one preprocessed object
[group: "Expert and maintenance"]
cluster-one input elbow_n="20" extra_dims="" resolutions="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--input "{{input}}" --elbow-n "{{elbow_n}}")
    if [ -n "{{extra_dims}}" ]; then
        args+=(--extra-dims "{{extra_dims}}")
    fi
    if [ -n "{{resolutions}}" ]; then
        args+=(--resolutions "{{resolutions}}")
    fi
    Rscript scripts/04-cluster.R "${args[@]}"

# Summarize full clustering grid
[group: "Expert and maintenance"]
summarize-clusters:
    Rscript scripts/05-summarize-clusters.R

# Summarize mg-selected clustering grid
[group: "Expert and maintenance"]
summarize-mg-selected elbow_n="20":
    Rscript scripts/08-summarize-mg-clusters.R --elbow-n {{elbow_n}}

# Plot full-dataset marker heatmap
[group: "Expert and maintenance"]
marker-heatmap dims="30" resolution="0.3" input="" layer="pflog" out_dir="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--dims "{{dims}}" --resolution "{{resolution}}" --layer "{{layer}}")
    if [ -n "{{input}}" ]; then
        args+=(--input "{{input}}")
    fi
    if [ -n "{{out_dir}}" ]; then
        args+=(--out-dir "{{out_dir}}")
    fi
    Rscript scripts/06-plot-marker-heatmap.R "${args[@]}"

# Plot cluster marker heatmaps
[group: "Expert and maintenance"]
cluster-marker-heatmaps dims="30" resolution="0.3" input="" layer="pflog" slot="data" n_perm="2000" out_dir="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--dims "{{dims}}" --resolution "{{resolution}}" --layer "{{layer}}" --slot "{{slot}}" --n-perm "{{n_perm}}")
    if [ -n "{{input}}" ]; then
        args+=(--input "{{input}}")
    fi
    if [ -n "{{out_dir}}" ]; then
        args+=(--out-dir "{{out_dir}}")
    fi
    Rscript scripts/10-plot-cluster-marker-heatmaps.R "${args[@]}"

# Run mg-selected marker ranking
[group: "Expert and maintenance"]
mg-markers input="" branch_tag="pflog_mg_selected_no_filter_cc" elbow_n="20" dims="20" resolution="0.5" assay="" layer="data" counts_layer="counts" top_n="5" min_pct="0.10" logfc_threshold="0.25" min_diff_pct="0" min_cells_group="3" cluster_map="" table_dir="" figure_dir="" confirm_no_merge="false" overwrite="false":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--branch-tag "{{branch_tag}}" --elbow-n "{{elbow_n}}" --dims "{{dims}}" --resolution "{{resolution}}" --layer "{{layer}}" --counts-layer "{{counts_layer}}" --top-n "{{top_n}}" --min-pct "{{min_pct}}" --logfc-threshold "{{logfc_threshold}}" --min-diff-pct "{{min_diff_pct}}" --min-cells-group "{{min_cells_group}}")
    if [ -n "{{input}}" ]; then args+=(--input "{{input}}"); fi
    if [ -n "{{assay}}" ]; then args+=(--assay "{{assay}}"); fi
    if [ -n "{{cluster_map}}" ]; then args+=(--cluster-map "{{cluster_map}}"); fi
    if [ -n "{{table_dir}}" ]; then args+=(--table-dir "{{table_dir}}"); fi
    if [ -n "{{figure_dir}}" ]; then args+=(--figure-dir "{{figure_dir}}"); fi
    if [ "{{confirm_no_merge}}" = "true" ]; then args+=(--confirm-no-merge); fi
    if [ "{{overwrite}}" = "true" ]; then args+=(--overwrite); fi
    Rscript scripts/11-find-mg-markers.R "${args[@]}"

# Run default mg-selected marker ranking without a cluster merge map
[group: "Expert and maintenance"]
mg-markers-no-merge:
    just mg-markers "" "pflog_mg_selected_no_filter_cc" "20" "20" "0.5" "" "data" "counts" "5" "0.10" "0.25" "0" "3" "" "" "" "true" "false"

# Plot mg-selected figures
[group: "Expert and maintenance"]
mg-figures input="" branch_tag="pflog_mg_selected_no_filter_cc" elbow_n="20" dims="20" resolution="0.5" layer="pflog" feature_list="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--branch-tag "{{branch_tag}}" --elbow-n "{{elbow_n}}" --dims "{{dims}}" --resolution "{{resolution}}" --layer "{{layer}}")
    if [ -n "{{input}}" ]; then args+=(--input "{{input}}"); fi
    if [ -n "{{feature_list}}" ]; then args+=(--feature-list "{{feature_list}}"); fi
    Rscript scripts/09-plot-mg-figures.R "${args[@]}"

# Run mg-selected DE and enrichment
[group: "Expert and maintenance"]
mg-de input="" cluster_column="cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5" condition_col="" control_label="" estim_label="" counts_layer="counts" deg_dir="" enrichment_dir="" lfc_shrink_type="apeglm" overwrite="false":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--cluster-column "{{cluster_column}}" --counts-layer "{{counts_layer}}" --lfc-shrink-type "{{lfc_shrink_type}}")
    if [ -n "{{input}}" ]; then args+=(--input "{{input}}"); fi
    if [ -n "{{condition_col}}" ]; then args+=(--condition-col "{{condition_col}}"); fi
    if [ -n "{{control_label}}" ]; then args+=(--control-label "{{control_label}}"); fi
    if [ -n "{{estim_label}}" ]; then args+=(--estim-label "{{estim_label}}"); fi
    if [ -n "{{deg_dir}}" ]; then args+=(--deg-dir "{{deg_dir}}"); fi
    if [ -n "{{enrichment_dir}}" ]; then args+=(--enrichment-dir "{{enrichment_dir}}"); fi
    if [ "{{overwrite}}" = "true" ]; then args+=(--overwrite); fi
    Rscript scripts/12-run-mg-de.R "${args[@]}"

# Re-run default mg-selected DE and replace existing outputs
[group: "Expert and maintenance"]
mg-de-overwrite:
    just mg-de "" "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5" "" "" "" "counts" "" "" "apeglm" "true"

# Render the Quarto notebook
[group: "Expert and maintenance"]
notebook:
    quarto render notebook/sc_analysis.qmd

# Run lightweight project tripwires
[group: "Expert and maintenance"]
tripwires:
    Rscript tools/run-tripwires.R
