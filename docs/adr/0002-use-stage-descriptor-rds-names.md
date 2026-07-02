# Use stage-descriptor names for saved Seurat objects

Pipeline code will use `sobj` as the working Seurat object name and save intermediate objects with `{step}_{descriptor}.rds` filenames. This keeps branch outputs sortable by analysis stage while preserving enough description to distinguish normalization branches, chosen clusterings, and downstream analysis checkpoints.
