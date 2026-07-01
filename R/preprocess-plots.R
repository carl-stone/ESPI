# Plots that go with scripts/preprocess-sobj.R

# TODO: Violin plot of nFeature_RNA, nCount_RNA, percent.mt, and percent.ribo for each sample.
# function name: splot_qc_metrics_violin (splot means "save plot"), indicating this is a plot byproduct function, not a returning function
# Make sure to explicitly group by sample because the sobj already has cluster identities that will
# be used by default with VlnPlot.
# 4 VlnPlots, one for each of the 4 metrics above
# Save as 2x2 grid of plots.
# Save as both PNG and PDF, 5 in by 5 in each.
# qc_metrics_violin.png/pdf

# TODO: HVG plot of the top n variable features
# function name: splot_hvg_scatter
# Args: sobj, n_top
# Saves scatter plot of mean expression vs variance for all genes, with the top n variable features highlighted in red.
# top_n <- head(VariableFeatures(sobj), n_top)
# p <- VariableFeaturePlot(sobj)
# p <- LabelPoints(plot = p, points = top_n, repel = TRUE)
# Save p as PNG and PDF, 4 in width by 3 in height
# hvg_scatter.png/pdf

# TODO: function to save VizDimLoadings plot of the top n PCs
# Must read reduction method (log1p vs pflog) and include that in the filename

# TODO: function to save ElbowPlot of the top n PCs
# Must read reduction method (log1p vs pflog) and include that in the filename

