# Exploratory plotting sandbox

`plot-sandbox.R` is an interactive workspace for expression, gene-pair,
cluster-focused, and sample-level detection exploration. It loads the configured
MG-selected object and current DESeq2 tables.

Run sections interactively from the repository root after
`devtools::load_all()`. It has no stable batch execution order and is not part
of `just run`. The script retains its existing Box output paths, including
`figures/random_pairs/`.

## Preserved prior exports

`outputs/all_gene_violins.pdf` and `outputs/all_ascl1_otx2_capb5_scatter.pdf`
were moved here from the repository root without regeneration. The sandbox
contains writers for these filenames, but the saved files differ in SHA-256
from the same-named current Box exports. Their exact producing run is not
established; retain them as prior exports, not current publication evidence.

Custom-marker manuscript work is separate: `scripts/custom-marker-plots.R`
writes to `notebook/figures/custom-markers/` (paths relative to the repo root).
