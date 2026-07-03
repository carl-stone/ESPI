cell_type_marker_genes <- list(
  muller_glia = c(
    "Rlbp1",
    "Glul",
    "Vim",
    "Slc1a3",
    "Sox9",
    "Hes1"
  ),
  activated_muller_glia = c(
    "Gfap",
    "Lcn2",
    "Serpina3n",
    "Ccn1",
    "Il6"
  ),
  proliferative = c(
    "Pcna",
    "Mcm2",
    "Mcm6",
    "Ccnd1",
    "Cdk4",
    "Cdk6"
  ),
  neurogenic_progenitor = c(
    "Ascl1",
    "Hes6",
    "Hes5",
    "Neurog2",
    "Dll1",
    "Dll3",
    "Neurod1",
    "Stmn1",
    "Stmn2"
  ),
  cone_bipolar = c(
    "Otx2",
    "Vsx2",
    "Cabp5",
    "Scgn",
    "Lhx4",
    "Grik1"
  ),
  rod_bipolar = c(
    "Prkca",
    "Grm6",
    "Car8"
  ),
  photoreceptor = c(
    "Crx",
    "Rcvrn",
    "Arr3",
    "Rho",
    "Esrrb",
    "Nrl",
    "Nr2e3",
    "Pde6g",
    "Rxrg"
  ),
  retinal_ganglion = c(
    "Chrna4",
    "Rbpms",
    "Thy1",
    "Sncg",
    "Slc17a6"
  ),
  microglia = c(
    "C1qa",
    "C1qb",
    "C1qc",
    "Cx3cr1",
    "Igf1"
  ),
  horizontal = c(
    "Onecut1",
    "Onecut2",
    "Lhx1",
    "Prox1",
    "Calb1",
    "Megf10",
    "Ntrk1"
  ),
  amacrine = c(
    "Tfap2a",
    "Tfap2b",
    "Pax6",
    "Calb2"
  )
)

cell_type_marker_labels <- c(
  muller_glia = "Müller glia",
  activated_muller_glia = "Activated Müller glia",
  proliferative = "Proliferative",
  neurogenic_progenitor = "Neurogenic progenitor",
  cone_bipolar = "Cone bipolar",
  rod_bipolar = "Rod bipolar",
  photoreceptor = "Photoreceptor",
  retinal_ganglion = "Retinal ganglion",
  microglia = "Microglia",
  horizontal = "Horizontal",
  amacrine = "Amacrine"
)

cell_type_marker_genes <- lapply(cell_type_marker_genes, unique)
if (
  length(setdiff(
    names(cell_type_marker_genes),
    names(cell_type_marker_labels)
  )) >
    0
) {
  stop(
    "Missing labels for cell types: ",
    paste(
      setdiff(names(cell_type_marker_genes), names(cell_type_marker_labels)),
      collapse = ", "
    )
  )
}

save(
  cell_type_marker_genes,
  file = "data/cell_type_marker_genes.rda",
  version = 2
)

save(
  cell_type_marker_labels,
  file = "data/cell_type_marker_labels.rda",
  version = 2
)
