human_cell_cycle_genes <- c(
  Seurat::cc.genes.updated.2019$s.genes,
  Seurat::cc.genes.updated.2019$g2m.genes
)

mouse <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "useast"
)

orthologs <- biomaRt::getBM(
  attributes = c(
    "external_gene_name",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type"
  ),
  filters = "with_hsapiens_homolog",
  values = TRUE,
  mart = mouse
)

mouse_cell_cycle_genes <- orthologs$external_gene_name[
  orthologs$hsapiens_homolog_associated_gene_name %in% human_cell_cycle_genes
]
mouse_cell_cycle_genes <- sort(unique(mouse_cell_cycle_genes))
mouse_cell_cycle_genes <- mouse_cell_cycle_genes[nzchar(mouse_cell_cycle_genes)]

save(mouse_cell_cycle_genes, file = "R/sysdata.rda", version = 2)
