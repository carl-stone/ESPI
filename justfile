set positional-arguments

# List available recipes
default:
    just --list

# Load package code
load:
    Rscript -e 'devtools::load_all(".", quiet = TRUE)'

# Update package documentation
document:
    Rscript -e 'devtools::document()'

# Rebuild README.md from README.Rmd
readme:
    Rscript -e 'devtools::build_readme()'

# Format R code with Air
format *paths:
    #!/usr/bin/env bash
    set -euo pipefail
    if [ "$#" -eq 0 ]; then
        air format R scripts
    else
        air format {{paths}}
    fi

# Run all preprocessing branches from legacy or counts-qc input
preprocess input_source="legacy" input="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--input-source "{{input_source}}")
    if [ -n "{{input}}" ]; then
        args=(--input "{{input}}")
    fi
    Rscript scripts/preprocess-all.R "${args[@]}"

# Run one preprocessing branch from legacy, counts-qc, or an explicit input
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
    Rscript scripts/preprocess-sobj.R "${args[@]}"

# Cluster every preprocessed object
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
    Rscript scripts/cluster-all.R "${args[@]}"

# Print cluster commands without running them
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
    Rscript scripts/cluster-all.R "${args[@]}"

# Cluster one preprocessed object
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
    Rscript scripts/cluster-sobj.R "${args[@]}"

# Summarize full clustering grid
summarize-clusters:
    Rscript scripts/summarize-cluster-grid.R

# Summarize mg-selected clustering grid
summarize-mg-selected elbow_n="20":
    Rscript scripts/summarize-mg-selected-grid.R --elbow-n {{elbow_n}}

# Plot full-dataset marker heatmap
marker-heatmap dims="50" resolution="0.3" input="" layer="pflog" out_dir="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--dims "{{dims}}" --resolution "{{resolution}}" --layer "{{layer}}")
    if [ -n "{{input}}" ]; then
        args+=(--input "{{input}}")
    fi
    if [ -n "{{out_dir}}" ]; then
        args+=(--out-dir "{{out_dir}}")
    fi
    Rscript scripts/big-heatmap-plot.R "${args[@]}"

# Plot cluster marker heatmaps
cluster-marker-heatmaps dims="50" resolution="0.3" input="" layer="pflog" slot="data" n_perm="2000" out_dir="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--dims "{{dims}}" --resolution "{{resolution}}" --layer "{{layer}}" --slot "{{slot}}" --n-perm "{{n_perm}}")
    if [ -n "{{input}}" ]; then
        args+=(--input "{{input}}")
    fi
    if [ -n "{{out_dir}}" ]; then
        args+=(--out-dir "{{out_dir}}")
    fi
    Rscript scripts/plot-cluster-marker-heatmaps.R "${args[@]}"

# Run mg-selected marker ranking
mg-markers input="" branch_tag="pflog_mg_selected_no_filter_cc" elbow_n="20" dims="30" resolution="0.3" assay="" layer="data" counts_layer="counts" top_n="5" min_pct="0.10" logfc_threshold="0.25" min_diff_pct="0" min_cells_group="3" cluster_map="" table_dir="" figure_dir="" confirm_no_merge="false" overwrite="false":
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
    Rscript scripts/find-markers-mg-selected.R "${args[@]}"

# Run default mg-selected marker ranking without a cluster merge map
mg-markers-no-merge:
    just mg-markers "" "pflog_mg_selected_no_filter_cc" "20" "30" "0.3" "" "data" "counts" "5" "0.10" "0.25" "0" "3" "" "" "" "true" "false"

# Plot mg-selected figures
mg-figures input="" branch_tag="pflog_mg_selected_no_filter_cc" elbow_n="20" dims="30" resolution="0.3" layer="pflog" feature_list="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--branch-tag "{{branch_tag}}" --elbow-n "{{elbow_n}}" --dims "{{dims}}" --resolution "{{resolution}}" --layer "{{layer}}")
    if [ -n "{{input}}" ]; then args+=(--input "{{input}}"); fi
    if [ -n "{{feature_list}}" ]; then args+=(--feature-list "{{feature_list}}"); fi
    Rscript scripts/plot-mg-selected-figures.R "${args[@]}"

# Run mg-selected DE and enrichment
mg-de input="" cluster_column="cluster_pflog_mg_selected_no_filter_cc_dims30_res0.3" condition_col="" control_label="" estim_label="" counts_layer="counts" deg_dir="" enrichment_dir="" lfc_shrink_type="normal" overwrite="false":
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
    Rscript scripts/run-mg-selected-de.R "${args[@]}"

# Re-run default mg-selected DE and replace existing outputs
mg-de-overwrite:
    just mg-de "" "cluster_pflog_mg_selected_no_filter_cc_dims30_res0.3" "" "" "" "counts" "" "" "normal" "true"

# Render the Quarto notebook
notebook:
    quarto render notebook/sc_analysis.qmd

# Run lightweight project tripwires
tripwires:
    Rscript tools/run-tripwires.R
